/*
 *
 *  This file is part of the SIRIUS library for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2023 Bright Giant GmbH
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with SIRIUS.
 *  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>
 */

package de.unijena.bioinf.storage.db.nosql.nitrite;

import com.fasterxml.jackson.core.Version;
import com.fasterxml.jackson.databind.JsonDeserializer;
import com.fasterxml.jackson.databind.JsonSerializer;
import com.fasterxml.jackson.databind.module.SimpleModule;
import com.google.common.collect.Iterables;
import de.unijena.bioinf.storage.db.nosql.Filter;
import de.unijena.bioinf.storage.db.nosql.Index;
import de.unijena.bioinf.storage.db.nosql.IndexType;
import de.unijena.bioinf.storage.db.nosql.*;
import org.apache.commons.lang3.tuple.Pair;
import org.dizitart.no2.*;
import org.dizitart.no2.filters.Filters;
import org.dizitart.no2.mapper.JacksonMapper;
import org.dizitart.no2.objects.ObjectFilter;
import org.dizitart.no2.objects.ObjectRepository;
import org.dizitart.no2.objects.filters.ObjectFilters;

import java.io.IOException;
import java.lang.reflect.Field;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantReadWriteLock;

public class NitriteDatabase implements Database<Document> {

    //todo I think this might cause even more problems. We should investigate this further.
    // But I think this cannot be executed on runtime on modern jdks. we need to take care of this as jvm parameters
    // Prevent illegal reflective access warnings
/*    static {
        if (!NitriteDatabase.class.getModule().isNamed()) {
            ClassLoader.class.getModule().addOpens(ClassLoader.class.getPackageName(), NitriteDatabase.class.getModule());
        }
    }*/
    protected Path file;

    // NITRITE
    private final Nitrite db;

    private final JacksonMapper nitriteMapper;

    private final Map<Class<?>, ObjectRepository<?>> repositories = Collections.synchronizedMap(new HashMap<>());

    private final Map<Class<?>, Set<String>> optionalRepoFields = Collections.synchronizedMap(new HashMap<>());

    private final Map<Class<?>, String> repoIdFields = Collections.synchronizedMap(new HashMap<>());

    private final Map<String, Set<String>> optionalCollectionFields = Collections.synchronizedMap(new HashMap<>());

    // LOCKS
    private final ReentrantReadWriteLock readWriteLock = new ReentrantReadWriteLock();
    private final ReentrantReadWriteLock.WriteLock writeLock = readWriteLock.writeLock();
    private final ReentrantReadWriteLock.ReadLock readLock = readWriteLock.readLock();
    private final ReentrantReadWriteLock stateLock = new ReentrantReadWriteLock();
    private final ReentrantReadWriteLock.WriteLock stateWriteLock = stateLock.writeLock();
    private final ReentrantReadWriteLock.ReadLock stateReadLock = stateLock.readLock();

    // STATE
    private boolean isClosed = false;

    public NitriteDatabase(Path file, Metadata meta) throws IOException {
        this.file = file;
        this.db = initDB(file, meta);
        this.initCollections(meta.collectionIndices);
        this.initRepositories(meta.repoIndices, meta.idFields);
        meta.optionalRepoFields.forEach((clazz, fields) -> optionalRepoFields.put(clazz, new HashSet<>(Arrays.asList(fields))));
        meta.optionalCollectionFields.forEach((collection, fields) -> optionalCollectionFields.put(collection, new HashSet<>(Arrays.asList(fields))));
        meta.idFields.forEach((clazz, pair) -> repoIdFields.put(clazz, pair.getKey()));
        try {

            Field cField = Nitrite.class.getDeclaredField("context");
            cField.setAccessible(true);
            NitriteContext context = (NitriteContext) cField.get(this.db);
            this.nitriteMapper = (JacksonMapper) context.getNitriteMapper();
        } catch (NoSuchFieldException | IllegalAccessException e) {
            throw new IOException(e);
        }
    }

    public JacksonMapper getJacksonMapper() {
        return this.nitriteMapper;
    }

    @SuppressWarnings("unchecked")
    private <T> void addSerializer(Metadata meta, SimpleModule module, Class<?> clazz, JsonSerializer<?> serializer) throws NoSuchFieldException {
        if (!meta.idFields.containsKey(clazz)) {
            Class<T> c = (Class<T>) clazz;
            JsonSerializer<T> s = (JsonSerializer<T>) serializer;
            module.addSerializer(c, s);
        }
    }

    @SuppressWarnings("unchecked")
    private <T> void addDeserializer(Metadata meta, SimpleModule module, Class<?> clazz, JsonDeserializer<?> deserializer) throws NoSuchFieldException {
        if (!meta.idFields.containsKey(clazz)) {
            Class<T> c = (Class<T>) clazz;
            JsonDeserializer<T> d = (JsonDeserializer<T>) deserializer;
            module.addDeserializer(c, d);
        }
    }

    @SuppressWarnings("unchecked")
    private <T> void addIdSerialization(Metadata meta, SimpleModule module, Class<T> clazz, String idField, boolean forceGenerateID) throws NoSuchFieldException, IllegalAccessException {
        if (!meta.serializers.containsKey(clazz)) {
            module.addSerializer(clazz, new NitriteIdMapperSerializer<>(clazz, idField, forceGenerateID, this));
        } else {
            module.addSerializer(clazz, new NitriteIdMapperSerializer<>(clazz, idField, forceGenerateID, (JsonSerializer<T>) meta.serializers.get(clazz)));
        }
        if (!meta.deserializers.containsKey(clazz)) {
            module.addDeserializer(clazz, new NitriteIdMapperDeserializer<>(clazz, idField, this));
        } else {
            module.addDeserializer(clazz, new NitriteIdMapperDeserializer<>(clazz, idField, (JsonDeserializer<T>) meta.deserializers.get(clazz)));
        }
    }

    private Nitrite initDB(Path file, Metadata meta) throws IOException {
        SimpleModule module = new SimpleModule("sirius-nitrite", Version.unknownVersion());
        try {
            for (Map.Entry<Class<?>, Pair<String, Boolean>> entry : meta.idFields.entrySet()) {
                addIdSerialization(meta, module, entry.getKey(), entry.getValue().getKey(), entry.getValue().getValue());
            }
            // FIXME conflict with id(de-)serializers!
            for (Map.Entry<Class<?>, JsonSerializer<?>> entry : meta.serializers.entrySet()) {
                addSerializer(meta, module, entry.getKey(), entry.getValue());
            }
            for (Map.Entry<Class<?>, JsonDeserializer<?>> entry : meta.deserializers.entrySet()) {
                addDeserializer(meta, module, entry.getKey(), entry.getValue());
            }
        } catch (NoSuchFieldException | IllegalAccessException e) {
            throw new IOException(e);
        }
        return Nitrite.builder().filePath(file.toFile()).registerModule(module).compressed().openOrCreate();
    }

    private void initCollections(Map<String, Index[]> collections) {
        for (String name : collections.keySet()) {
            NitriteCollection collection = this.db.getCollection(name);
            initIndex(collections.get(name), collection);
        }
    }

    private void initRepositories(Map<Class<?>, Index[]> repositories, Map<Class<?>, Pair<String, Boolean>> idFields) throws IOException {
        for (Class<?> clazz : repositories.keySet()) {
            ObjectRepository<?> repository = this.db.getRepository(clazz);
            this.repositories.put(clazz, repository);
            initIndex(repositories.get(clazz), repository);
            if (idFields.containsKey(clazz)) {
                initIdField(clazz, idFields.get(clazz).getKey(), repository);
            }
        }
    }

    private void initIdField(Class<?> clazz, String idFieldName, ObjectRepository<?> repository) throws IOException {
        try {
            Field idField = clazz.getDeclaredField(idFieldName);
            Field fieldField = repository.getClass().getDeclaredField("idField");
            fieldField.setAccessible(true);
            fieldField.set(repository, idField);
            if (!repository.hasIndex(idFieldName)) {
                repository.createIndex(idFieldName, IndexOptions.indexOptions(org.dizitart.no2.IndexType.Unique));
            }
        } catch (NoSuchFieldException | IllegalAccessException e) {
            throw new IOException(e);
        }
    }

    private <Repo extends PersistentCollection<?>> void initIndex(Index[] indices, Repo repository) {
        Collection<org.dizitart.no2.Index> repoIndices = repository.listIndices();
        List<org.dizitart.no2.Index> toDrop = new ArrayList<>();
        List<Index> toBuild = new ArrayList<>();
        for (Index index : indices) {
            found:
            {
                for (org.dizitart.no2.Index repoIndex : repoIndices) {
                    if (Objects.equals(repoIndex.getField(), index.getField())) {
                        org.dizitart.no2.IndexType repoType = repoIndex.getIndexType();
                        IndexType iType = index.getType();
                        if ((repoType == org.dizitart.no2.IndexType.Unique && iType != IndexType.UNIQUE) ||
                                (repoType == org.dizitart.no2.IndexType.NonUnique && iType != IndexType.NON_UNIQUE) ||
                                (repoType == org.dizitart.no2.IndexType.Fulltext && iType != IndexType.FULL_TEXT)) {
                            toDrop.add(repoIndex);
                            toBuild.add(index);
                        }
                        break found;
                    }
                }
                toBuild.add(index);
            }
        }

        for (org.dizitart.no2.Index index : toDrop) {
            repository.dropIndex(index.getField());
        }

        for (Index index : toBuild) {
            switch (index.getType()) {
                case UNIQUE:
                    repository.createIndex(index.getField(), IndexOptions.indexOptions(org.dizitart.no2.IndexType.Unique));
                    break;
                case NON_UNIQUE:
                    repository.createIndex(index.getField(), IndexOptions.indexOptions(org.dizitart.no2.IndexType.NonUnique));
                    break;
                case FULL_TEXT:
                    repository.createIndex(index.getField(), IndexOptions.indexOptions(org.dizitart.no2.IndexType.Fulltext));
                    break;
            }
        }
    }

    @Override
    public void close() {
        stateWriteLock.lock();
        try {
            this.isClosed = true;
            this.db.close();
        } finally {
            stateWriteLock.unlock();
        }
    }

    private <T> T callIfOpen(Callable<T> callable) throws IOException {
        stateReadLock.lock();
        if (this.isClosed) {
            throw new IOException("Nitrite database is closed!");
        }
        try {
            return callable.call();
        } catch (IOException e) {
            throw e;
        } catch (Exception e) {
            throw new RuntimeException(e);
        } finally {
            stateReadLock.unlock();
        }
    }

    private <T> T read(Callable<T> callable) throws IOException {
        return this.callIfOpen(() -> {
            readLock.lock();
            try {
                return callable.call();
            } finally {
                readLock.unlock();
            }
        });
    }

    private <T> T write(Callable<T> callable) throws IOException {
        return this.callIfOpen(() -> {
            writeLock.lock();
            try {
                return callable.call();
            } finally {
                writeLock.unlock();
            }
        });
    }

    @SuppressWarnings("unchecked")
    private <T> ObjectRepository<T> getRepository(Class<T> clazz) throws IOException {
        if (!this.repositories.containsKey(clazz)) {
            throw new IOException(clazz + " is not registered.");
        }
        return (ObjectRepository<T>) this.repositories.get(clazz);
    }

    @SuppressWarnings("unchecked")
    private <T> ObjectRepository<T> getRepository(T object) throws IOException {
        if (!this.repositories.containsKey(object.getClass())) {
            throw new IOException(object.getClass() + " is not registered.");
        }
        return (ObjectRepository<T>) this.repositories.get(object.getClass());
    }

    @SuppressWarnings("unchecked")
    private <T> Pair<T[], ObjectRepository<T>> getRepository(Iterable<T> objects) throws IOException {
        Collection<T> collection = new ArrayList<>();
        Iterables.addAll(collection, objects);
        if (collection.isEmpty()) {
            return null;
        }
        T[] arr = (T[]) collection.toArray();
        Class<T> clazz = (Class<T>) arr[0].getClass();
        if (!this.repositories.containsKey(clazz)) {
            throw new IOException(clazz + " is not registered.");
        }
        return Pair.of(arr, (ObjectRepository<T>) this.repositories.get(clazz));
    }

    private NitriteCollection getCollection(String name) throws IOException {
        //collection of repos are created on demand, so we have to check both
        if (!db.hasCollection(name) && !db.listRepositories().contains(name)) {
            throw new IOException(name + " is not registered.");
        }

        return db.getCollection(name);
    }

    private <T> Iterable<T> maybeProject(Class<T> clazz, org.dizitart.no2.objects.Cursor<T> cursor, String[] withOptionalFields) throws IOException {
        if (optionalRepoFields.containsKey(clazz)) {
            Set<String> omittedFields = new HashSet<>(optionalRepoFields.get(clazz));
            omittedFields.removeAll(new HashSet<>(Arrays.asList(withOptionalFields)));
            if (omittedFields.size() > 0) {
                return new ProjectingIterable<>(clazz, cursor, omittedFields, nitriteMapper);
            }
        }
        return cursor;
    }

    private Iterable<Document> maybeProjectDocuments(String collectionName, Cursor cursor, String[] withOptionalFields) {
        if (optionalCollectionFields.containsKey(collectionName)) {
            Set<String> omittedFields = new HashSet<>(optionalCollectionFields.get(collectionName));
            omittedFields.removeAll(new HashSet<>(Arrays.asList(withOptionalFields)));
            if (omittedFields.size() > 0) {
                return new ProjectingDocumentIterable(cursor, omittedFields);
            }
        }
        return cursor;
    }

    @Override
    public Path location() {
        return file;
    }

    @Override
    @SuppressWarnings("unchecked")
    public <T> int insert(T object) throws IOException {
        return this.write(() -> {
            ObjectRepository<T> repo = this.getRepository(object);
            return repo.insert(object).getAffectedCount();
        });
    }

    @Override
    public <T> int insertAll(Iterable<T> objects) throws IOException {
        return this.write(() -> {
            Pair<T[], ObjectRepository<T>> pair = this.getRepository(objects);
            if (pair == null) {
                return 0;
            }
            return pair.getRight().insert(pair.getLeft()).getAffectedCount();
        });
    }

    @Override
    public int insert(String collectionName, Document document) throws IOException {
        return this.write(() -> {
            NitriteCollection collection = this.getCollection(collectionName);
            return collection.insert(document).getAffectedCount();
        });
    }

    @Override
    public int insertAll(String collectionName, Iterable<Document> documents) throws IOException {
        return this.write(() -> {
            NitriteCollection collection = this.getCollection(collectionName);
            Document[] docs = Iterables.toArray(documents, Document.class);
            return collection.insert(docs).getAffectedCount();
        });
    }

    @Override
    public <T> int upsert(T object) throws IOException {
        return this.write(() -> {
            ObjectRepository<T> repo = this.getRepository(object);
            return repo.update(object, true).getAffectedCount();
        });
    }

    @Override
    public <T> int upsertAll(Iterable<T> objects) throws IOException {
        return this.write(() -> {
            Pair<T[], ObjectRepository<T>> pair = this.getRepository(objects);
            if (pair == null) {
                return 0;
            }
            int count = 0;
            for (T o : pair.getLeft()) {
                count += pair.getRight().update(o, true).getAffectedCount();
            }
            return count;
        });
    }

    @Override
    public int upsert(String collectionName, Document document) throws IOException {
        return this.write(() -> {
            NitriteCollection collection = this.getCollection(collectionName);
            return collection.update(document, true).getAffectedCount();
        });
    }

    @Override
    public int upsertAll(String collectionName, Iterable<Document> documents) throws IOException {
        return this.write(() -> {
            NitriteCollection collection = this.getCollection(collectionName);
            int count = 0;
            for (Document doc : documents) {
                count += collection.update(doc, true).getAffectedCount();
            }
            return count;
        });
    }

    @Override
    public <T> T getById(long id, Class<T> clazz, String... withOptionalFields) throws IOException {
        return this.read(() -> {
            ObjectRepository<T> repo = this.getRepository(clazz);
            if (optionalRepoFields.containsKey(clazz)) {
                Set<String> omittedFields = new HashSet<>(optionalRepoFields.get(clazz));
                omittedFields.removeAll(new HashSet<>(Arrays.asList(withOptionalFields)));
                if (omittedFields.size() > 0) {
                    Field cField = repo.getClass().getDeclaredField("collection");
                    cField.setAccessible(true);
                    NitriteCollection collection = (NitriteCollection) cField.get(repo);
                    Document document = collection.getById(NitriteId.createId(id));
                    return nitriteMapper.asObject(ProjectingDocumentIterable.project(document, omittedFields), clazz);
                } else {
                    return repo.getById(NitriteId.createId(id));
                }
            } else {
                return repo.getById(NitriteId.createId(id));
            }
        });
    }

    @Override
    public Document getById(String collectionName, long id, String... withOptionalFields) throws IOException {
        return this.read(() -> {
            NitriteCollection collection = this.getCollection(collectionName);
            Document document = collection.getById(NitriteId.createId(id));
            if (optionalCollectionFields.containsKey(collectionName)) {
                Set<String> omittedFields = new HashSet<>(optionalCollectionFields.get(collectionName));
                omittedFields.removeAll(new HashSet<>(Arrays.asList(withOptionalFields)));
                if (omittedFields.size() > 0) {
                    return ProjectingDocumentIterable.project(document, omittedFields);
                } else {
                    return document;
                }
            } else {
                return document;
            }
        });
    }

    @Override
    public <T> Iterable<T> find(Filter filter, Class<T> clazz, String... withOptionalFields) throws IOException {
        return this.read(() -> {
            ObjectRepository<T> repo = this.getRepository(clazz);
            ObjectFilter of = getObjectFilter(filter);
            return maybeProject(clazz, repo.find(of), withOptionalFields);
        });
    }

    @Override
    public <T> Iterable<T> find(Filter filter, Class<T> clazz, int offset, int pageSize, String... withOptionalFields) throws IOException {
        return this.read(() -> {
            ObjectRepository<T> repo = this.getRepository(clazz);
            ObjectFilter of = getObjectFilter(filter);
            return maybeProject(clazz, repo.find(of, FindOptions.limit(offset, pageSize)), withOptionalFields);
        });
    }

    @Override
    public <T> Iterable<T> find(Filter filter, Class<T> clazz, String sortField, SortOrder sortOrder, String... withOptionalFields) throws IOException {
        return this.read(() -> {
            ObjectRepository<T> repo = this.getRepository(clazz);
            ObjectFilter of = getObjectFilter(filter);
            return maybeProject(clazz, repo.find(of, FindOptions.sort(sortField, (sortOrder == SortOrder.ASCENDING) ? org.dizitart.no2.SortOrder.Ascending : org.dizitart.no2.SortOrder.Descending)), withOptionalFields);
        });
    }

    @Override
    public <T> Iterable<T> find(Filter filter, Class<T> clazz, int offset, int pageSize, String sortField, SortOrder sortOrder, String... withOptionalFields) throws IOException {
        return this.read(() -> {
            ObjectRepository<T> repo = this.getRepository(clazz);
            ObjectFilter of = getObjectFilter(filter);
            return maybeProject(clazz, repo.find(of, FindOptions.sort(sortField, (sortOrder == SortOrder.ASCENDING) ? org.dizitart.no2.SortOrder.Ascending : org.dizitart.no2.SortOrder.Descending).thenLimit(offset, pageSize)), withOptionalFields);
        });
    }

    @Override
    public Iterable<Document> find(String collectionName, Filter filter, String... withOptionalFields) throws IOException {
        return this.read(() -> {
            NitriteCollection collection = this.getCollection(collectionName);
            org.dizitart.no2.Filter f = getFilter(filter);
            return maybeProjectDocuments(collectionName, collection.find(f), withOptionalFields);
        });
    }

    @Override
    public Iterable<Document> find(String collectionName, Filter filter, int offset, int pageSize, String... withOptionalFields) throws IOException {
        return this.read(() -> {
            NitriteCollection collection = this.getCollection(collectionName);
            org.dizitart.no2.Filter f = getFilter(filter);
            return maybeProjectDocuments(collectionName, collection.find(f, FindOptions.limit(offset, pageSize)), withOptionalFields);
        });
    }

    @Override
    public Iterable<Document> find(String collectionName, Filter filter, String sortField, SortOrder sortOrder, String... withOptionalFields) throws IOException {
        return this.read(() -> {
            NitriteCollection collection = this.getCollection(collectionName);
            org.dizitart.no2.Filter f = getFilter(filter);
            return maybeProjectDocuments(collectionName, collection.find(f, FindOptions.sort(sortField, (sortOrder == SortOrder.ASCENDING) ? org.dizitart.no2.SortOrder.Ascending : org.dizitart.no2.SortOrder.Descending)), withOptionalFields);
        });
    }

    @Override
    public Iterable<Document> find(String collectionName, Filter filter, int offset, int pageSize, String sortField, SortOrder sortOrder, String... withOptionalFields) throws IOException {
        return this.read(() -> {
            NitriteCollection collection = this.getCollection(collectionName);
            org.dizitart.no2.Filter f = getFilter(filter);
            return maybeProjectDocuments(collectionName, collection.find(f, FindOptions.sort(sortField, (sortOrder == SortOrder.ASCENDING) ? org.dizitart.no2.SortOrder.Ascending : org.dizitart.no2.SortOrder.Descending).thenLimit(offset, pageSize)), withOptionalFields);
        });
    }

    @Override
    public <T> Iterable<T> findAll(Class<T> clazz, String... withOptionalFields) throws IOException {
        return this.read(() -> {
            ObjectRepository<T> repo = this.getRepository(clazz);
            return maybeProject(clazz, repo.find(), withOptionalFields);
        });
    }

    @Override
    public <T> Iterable<T> findAll(Class<T> clazz, int offset, int pageSize, String... withOptionalFields) throws IOException {
        return this.read(() -> {
            ObjectRepository<T> repo = this.getRepository(clazz);
            return maybeProject(clazz, repo.find(FindOptions.limit(offset, pageSize)), withOptionalFields);
        });
    }

    @Override
    public <T> Iterable<T> findAll(Class<T> clazz, String sortField, SortOrder sortOrder, String... withOptionalFields) throws IOException {
        return this.read(() -> {
            ObjectRepository<T> repo = this.getRepository(clazz);
            return maybeProject(clazz, repo.find(
                    FindOptions.sort(sortField, (sortOrder == SortOrder.ASCENDING) ? org.dizitart.no2.SortOrder.Ascending : org.dizitart.no2.SortOrder.Descending)
            ), withOptionalFields);
        });
    }

    @Override
    public <T> Iterable<T> findAll(Class<T> clazz, int offset, int pageSize, String sortField, SortOrder sortOrder, String... withOptionalFields) throws IOException {
        return this.read(() -> {
            ObjectRepository<T> repo = this.getRepository(clazz);
            return maybeProject(clazz, repo.find(
                    FindOptions.sort(sortField, (sortOrder == SortOrder.ASCENDING) ? org.dizitart.no2.SortOrder.Ascending : org.dizitart.no2.SortOrder.Descending).thenLimit(offset, pageSize)
            ), withOptionalFields);
        });
    }

    @Override
    public Iterable<Document> findAll(String collectionName, String... withOptionalFields) throws IOException {
        return this.read(() -> {
            NitriteCollection collection = this.getCollection(collectionName);
            return maybeProjectDocuments(collectionName, collection.find(), withOptionalFields);
        });
    }

    @Override
    public Iterable<Document> findAll(String collectionName, int offset, int pageSize, String... withOptionalFields) throws IOException {
        return this.read(() -> {
            NitriteCollection collection = this.getCollection(collectionName);
            return maybeProjectDocuments(collectionName, collection.find(FindOptions.limit(offset, pageSize)), withOptionalFields);
        });
    }

    @Override
    public Iterable<Document> findAll(String collectionName, String sortField, SortOrder sortOrder, String... withOptionalFields) throws IOException {
        return this.read(() -> {
            NitriteCollection collection = this.getCollection(collectionName);
            return maybeProjectDocuments(collectionName, collection.find(
                    FindOptions.sort(sortField, (sortOrder == SortOrder.ASCENDING) ? org.dizitart.no2.SortOrder.Ascending : org.dizitart.no2.SortOrder.Descending)
            ), withOptionalFields);
        });
    }

    @Override
    public Iterable<Document> findAll(String collectionName, int offset, int pageSize, String sortField, SortOrder sortOrder, String... withOptionalFields) throws IOException {
        return this.read(() -> {
            NitriteCollection collection = this.getCollection(collectionName);
            return maybeProjectDocuments(collectionName, collection.find(
                    FindOptions.sort(sortField, (sortOrder == SortOrder.ASCENDING) ? org.dizitart.no2.SortOrder.Ascending : org.dizitart.no2.SortOrder.Descending).thenLimit(offset, pageSize)
            ), withOptionalFields);
        });
    }

    @Override
    public <T> T injectOptionalFields(T object, String... optionalFields) throws IOException {
        NitriteCollection collection;
        try {
            ObjectRepository<T> repository = getRepository(object);
            Field cField = repository.getClass().getDeclaredField("collection");
            cField.setAccessible(true);
            collection = (NitriteCollection) cField.get(repository);
        } catch (NoSuchFieldException | IllegalAccessException e) {
            throw new IOException(e);
        }
        if (!repoIdFields.containsKey(object.getClass())) {
            throw new IOException("Object has no ID field.");
        }
        return InjectingIterable.inject(object, new HashSet<>(Arrays.asList(optionalFields)), collection, repoIdFields.get(object.getClass()), nitriteMapper);
    }

    @Override
    public Document injectOptionalFields(String collectionName, Document document, String... optionalFields) throws IOException {
        return InjectingDocumentIterable.inject(document, new HashSet<>(Arrays.asList(optionalFields)), getCollection(collectionName));
    }

    @Override
    public <T> Iterable<T> injectOptionalFields(Class<T> clazz, Iterable<T> objects, String... optionalFields) throws IOException {
        if (objects instanceof org.dizitart.no2.objects.Cursor<T>) {
            return maybeProject(clazz, (org.dizitart.no2.objects.Cursor<T>) objects, optionalFields);
        } else if (objects instanceof ProjectingIterable<T>) {
            ((ProjectingIterable<T>) objects).withOptionalFields(new HashSet<>(Arrays.asList(optionalFields)));
            return objects;
        } else {
            if (!repoIdFields.containsKey(clazz)) {
                throw new IOException("Object has no ID field.");
            }
            return new InjectingIterable<>(objects, new HashSet<>(Arrays.asList(optionalFields)), getRepository(clazz), repoIdFields.get(clazz), nitriteMapper);
        }
    }

    @Override
    public Iterable<Document> injectOptionalFields(String collectionName, Iterable<Document> documents, String... optionalFields) throws IOException {
        if (documents instanceof Cursor) {
            return maybeProjectDocuments(collectionName, (Cursor) documents, optionalFields);
        } else if (documents instanceof ProjectingDocumentIterable) {
            ((ProjectingDocumentIterable) documents).withOptionalFields(new HashSet<>(Arrays.asList(optionalFields)));
            return documents;
        } else {
            return new InjectingDocumentIterable(documents, new HashSet<>(Arrays.asList(optionalFields)), getCollection(collectionName));
        }
    }

    @Override
    public <P, C> Iterable<P> joinAllChildren(Class<C> childClass, Iterable<P> parents, String localField, String foreignField, String targetField, String... withOptionalChildFields) throws IOException {
        return new JoinedReflectionIterable<>(
                childClass,
                parents,
                (localObject) -> {
                    try {
                        org.dizitart.no2.objects.Cursor<C> objectCursor = (org.dizitart.no2.objects.Cursor<C>) find(new Filter().eq(foreignField, localObject), childClass);
                        Field cField = objectCursor.getClass().getDeclaredField("cursor");
                        cField.setAccessible(true);
                        return maybeProjectDocuments(childClass.getName(), (Cursor) cField.get(objectCursor), withOptionalChildFields);
                    } catch (IOException | IllegalAccessException | NoSuchFieldException e) {
                        throw new RuntimeException(e);
                    }
                },
                localField,
                targetField,
                nitriteMapper
        );
    }

    @Override
    public <P, C> Iterable<P> joinChildren(Class<C> childClass, Filter childFilter, Iterable<P> parents, String localField, String foreignField, String targetField, String... withOptionalChildFields) throws IOException {
        return new JoinedReflectionIterable<>(
                childClass,
                parents,
                (localObject) -> {
                    try {
                        Filter cFilter = new Filter().and().eq(foreignField, localObject);
                        cFilter.filterChain.addAll(childFilter.filterChain);
                        org.dizitart.no2.objects.Cursor<C> objectCursor = (org.dizitart.no2.objects.Cursor<C>) find(cFilter, childClass);
                        Field cField = objectCursor.getClass().getDeclaredField("cursor");
                        cField.setAccessible(true);
                        return maybeProjectDocuments(childClass.getName(), (Cursor) cField.get(objectCursor), withOptionalChildFields);
                    } catch (IOException | IllegalAccessException | NoSuchFieldException e) {
                        throw new RuntimeException(e);
                    }
                },
                localField,
                targetField,
                nitriteMapper
        );
    }

    @Override
    public <T, P, C> Iterable<T> joinAllChildren(Class<T> targetClass, Class<C> childClass, Iterable<P> parents, String localField, String foreignField, String targetField, String... withOptionalChildFields) throws IOException {
        return new JoinedIterable<>(
                targetClass,
                parents,
                (localObject) -> {
                    try {
                        org.dizitart.no2.objects.Cursor<C> objectCursor = (org.dizitart.no2.objects.Cursor<C>) find(new Filter().eq(foreignField, localObject), childClass);
                        Field cField = objectCursor.getClass().getDeclaredField("cursor");
                        cField.setAccessible(true);
                        return maybeProjectDocuments(childClass.getName(), (Cursor) cField.get(objectCursor), withOptionalChildFields);
                    } catch (IOException | IllegalAccessException | NoSuchFieldException e) {
                        throw new RuntimeException(e);
                    }
                },
                localField,
                targetField,
                nitriteMapper
        );
    }

    @Override
    public <T, P, C> Iterable<T> joinChildren(Class<T> targetClass, Class<C> childClass, Filter childFilter, Iterable<P> parents, String localField, String foreignField, String targetField, String... withOptionalChildFields) throws IOException {
        return new JoinedIterable<>(
                targetClass,
                parents,
                (localObject) -> {
                    try {
                        Filter cFilter = new Filter().and().eq(foreignField, localObject);
                        cFilter.filterChain.addAll(childFilter.filterChain);
                        org.dizitart.no2.objects.Cursor<C> objectCursor = (org.dizitart.no2.objects.Cursor<C>) find(cFilter, childClass);
                        Field cField = objectCursor.getClass().getDeclaredField("cursor");
                        cField.setAccessible(true);
                        return maybeProjectDocuments(childClass.getName(), (Cursor) cField.get(objectCursor), withOptionalChildFields);
                    } catch (IOException | IllegalAccessException | NoSuchFieldException e) {
                        throw new RuntimeException(e);
                    }
                },
                localField,
                targetField,
                nitriteMapper
        );
    }

    @Override
    public Iterable<Document> joinAllChildren(String childCollectionName, Iterable<Document> parents, String localField, String foreignField, String targetField, String... withOptionalChildFields) throws IOException {
        return () -> new JoinedDocumentIterator(
                parents,
                (localObject) -> {
                    try {
                        return find(childCollectionName, new Filter().eq(foreignField, localObject), withOptionalChildFields);
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                },
                localField,
                targetField
        );
    }

    @Override
    public Iterable<Document> joinChildren(String childCollectionName, Filter childFilter, Iterable<Document> parents, String localField, String foreignField, String targetField, String... withOptionalChildFields) throws IOException {
        return () -> new JoinedDocumentIterator(
                parents,
                (localObject) -> {
                    try {
                        Filter cFilter = new Filter().and().eq(foreignField, localObject);
                        cFilter.filterChain.addAll(childFilter.filterChain);
                        return find(childCollectionName, cFilter, withOptionalChildFields);
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                },
                localField,
                targetField
        );
    }

    @Override
    public <T> int count(Filter filter, Class<T> clazz) throws IOException {
        return this.read(() -> {
            ObjectRepository<T> repo = this.getRepository(clazz);
            ObjectFilter of = getObjectFilter(filter);
            return repo.find(of).size();
        });
    }

    @Override
    public <T> int count(Filter filter, Class<T> clazz, int offset, int pageSize) throws IOException {
        return this.read(() -> {
            ObjectRepository<T> repo = this.getRepository(clazz);
            ObjectFilter of = getObjectFilter(filter);
            return repo.find(of, FindOptions.limit(offset, pageSize)).size();
        });
    }

    @Override
    public <T> int countAll(Class<T> clazz) throws IOException {
        return this.read(() -> {
            ObjectRepository<T> repo = this.getRepository(clazz);
            return repo.find().totalCount();
        });
    }

    @Override
    public int count(String collectionName, Filter filter) throws IOException {
        return this.read(() -> {
            NitriteCollection collection = this.getCollection(collectionName);
            org.dizitart.no2.Filter f = getFilter(filter);
            return collection.find(f).size();
        });
    }

    @Override
    public int count(String collectionName, Filter filter, int offset, int pageSize) throws IOException {
        return this.read(() -> {
            NitriteCollection collection = this.getCollection(collectionName);
            org.dizitart.no2.Filter f = getFilter(filter);
            return collection.find(f, FindOptions.limit(offset, pageSize)).size();
        });
    }

    @Override
    public int countAll(String collectionName) throws IOException {
        return this.read(() -> {
            NitriteCollection collection = this.getCollection(collectionName);
            return collection.find().totalCount();
        });
    }

    @Override
    public <T> int remove(T object) throws IOException {
        return this.write(() -> {
            ObjectRepository<T> repo = this.getRepository(object);
            return repo.remove(object).getAffectedCount();
        });
    }

    @Override
    public <T> int removeAll(Iterable<T> objects) throws IOException {
        return this.write(() -> {
            Pair<T[], ObjectRepository<T>> pair = this.getRepository(objects);
            if (pair == null) {
                return 0;
            }
            int count = 0;
            for (T o : pair.getLeft()) {
                count += pair.getRight().remove(o).getAffectedCount();
            }
            return count;
        });
    }

    @Override
    public <T> int removeAll(Filter filter, Class<T> clazz) throws IOException {
        return this.write(() -> {
            ObjectRepository<T> repo = this.getRepository(clazz);
            ObjectFilter of = getObjectFilter(filter);
            return repo.remove(of).getAffectedCount();
        });
    }

    @Override
    public int remove(String collectionName, Document document) throws IOException {
        return this.write(() -> {
            NitriteCollection collection = this.getCollection(collectionName);
            return collection.remove(document).getAffectedCount();
        });
    }

    @Override
    public int removeAll(String collectionName, Iterable<Document> documents) throws IOException {
        return this.write(() -> {
            NitriteCollection collection = this.getCollection(collectionName);
            int count = 0;
            for (Document doc : documents) {
                count += collection.remove(doc).getAffectedCount();
            }
            return count;
        });
    }

    @Override
    public int removeAll(String collectionName, Filter filter) throws IOException {
        return this.write(() -> {
            NitriteCollection collection = this.getCollection(collectionName);
            org.dizitart.no2.Filter f = getFilter(filter);
            return collection.remove(f).getAffectedCount();
        });
    }

    private ObjectFilter[] getAllObjectFilterChildren(Deque<Filter.FilterElement> filterChain) {
        List<ObjectFilter> children = new ArrayList<>();
        while (!filterChain.isEmpty()) {
            ObjectFilter next = getObjectFilter(filterChain);
            if (next == ObjectFilters.ALL) {
                break;
            }
            children.add(next);
        }
        ObjectFilter[] arr = new ObjectFilter[children.size()];
        return children.toArray(arr);
    }

    private ObjectFilter getObjectFilter(Deque<Filter.FilterElement> filterChain) {
        if (filterChain.isEmpty()) {
            return ObjectFilters.ALL;
        }
        Filter.FilterElement element = filterChain.pop();
        switch (element.filterType) {
            case AND:
                return ObjectFilters.and(getAllObjectFilterChildren(filterChain));
            case OR:
                return ObjectFilters.or(getAllObjectFilterChildren(filterChain));
            case NOT:
                return ObjectFilters.not(getObjectFilter(filterChain));
            case EQ:
                return ObjectFilters.eq(
                        ((Filter.FieldFilterElement) element).field,
                        ((Filter.FieldFilterElement) element).values[0]
                );
            case GT:
                return ObjectFilters.gt(
                        ((Filter.FieldFilterElement) element).field,
                        ((Filter.FieldFilterElement) element).values[0]
                );
            case GTE:
                return ObjectFilters.gte(
                        ((Filter.FieldFilterElement) element).field,
                        ((Filter.FieldFilterElement) element).values[0]
                );
            case LT:
                return ObjectFilters.lt(
                        ((Filter.FieldFilterElement) element).field,
                        ((Filter.FieldFilterElement) element).values[0]
                );
            case LTE:
                return ObjectFilters.lte(
                        ((Filter.FieldFilterElement) element).field,
                        ((Filter.FieldFilterElement) element).values[0]
                );
            case TEXT:
                return ObjectFilters.text(
                        ((Filter.FieldFilterElement) element).field,
                        (String) ((Filter.FieldFilterElement) element).values[0]
                );
            case REGEX:
                return ObjectFilters.regex(
                        ((Filter.FieldFilterElement) element).field,
                        (String) ((Filter.FieldFilterElement) element).values[0]
                );
            case IN:
                return ObjectFilters.in(
                        ((Filter.FieldFilterElement) element).field,
                        ((Filter.FieldFilterElement) element).values
                );
            case NOT_IN:
                return ObjectFilters.notIn(
                        ((Filter.FieldFilterElement) element).field,
                        ((Filter.FieldFilterElement) element).values
                );
            case ELEM_MATCH:
                return ObjectFilters.elemMatch(
                        ((Filter.FieldFilterElement) element).field,
                        getObjectFilter(filterChain)
                );
        }
        return ObjectFilters.ALL;
    }

    private org.dizitart.no2.Filter[] getAllFilterChildren(Deque<Filter.FilterElement> filterChain) {
        List<org.dizitart.no2.Filter> children = new ArrayList<>();
        while (!filterChain.isEmpty()) {
            org.dizitart.no2.Filter next = getFilter(filterChain);
            if (next == Filters.ALL) {
                break;
            }
            children.add(next);
        }
        org.dizitart.no2.Filter[] arr = new org.dizitart.no2.Filter[children.size()];
        return children.toArray(arr);
    }

    private org.dizitart.no2.Filter getFilter(Deque<Filter.FilterElement> filterChain) {
        if (filterChain.isEmpty()) {
            return Filters.ALL;
        }
        Filter.FilterElement element = filterChain.pop();
        switch (element.filterType) {
            case AND:
                return Filters.and(getAllFilterChildren(filterChain));
            case OR:
                return Filters.or(getAllFilterChildren(filterChain));
            case NOT:
                return Filters.not(getFilter(filterChain));
            case EQ:
                return Filters.eq(
                        ((Filter.FieldFilterElement) element).field,
                        ((Filter.FieldFilterElement) element).values[0]
                );
            case GT:
                return Filters.gt(
                        ((Filter.FieldFilterElement) element).field,
                        ((Filter.FieldFilterElement) element).values[0]
                );
            case GTE:
                return Filters.gte(
                        ((Filter.FieldFilterElement) element).field,
                        ((Filter.FieldFilterElement) element).values[0]
                );
            case LT:
                return Filters.lt(
                        ((Filter.FieldFilterElement) element).field,
                        ((Filter.FieldFilterElement) element).values[0]
                );
            case LTE:
                return Filters.lte(
                        ((Filter.FieldFilterElement) element).field,
                        ((Filter.FieldFilterElement) element).values[0]
                );
            case TEXT:
                return Filters.text(
                        ((Filter.FieldFilterElement) element).field,
                        (String) ((Filter.FieldFilterElement) element).values[0]
                );
            case REGEX:
                return Filters.regex(
                        ((Filter.FieldFilterElement) element).field,
                        (String) ((Filter.FieldFilterElement) element).values[0]
                );
            case IN:
                return Filters.in(
                        ((Filter.FieldFilterElement) element).field,
                        ((Filter.FieldFilterElement) element).values
                );
            case NOT_IN:
                return Filters.notIn(
                        ((Filter.FieldFilterElement) element).field,
                        ((Filter.FieldFilterElement) element).values
                );
            case ELEM_MATCH:
                return Filters.elemMatch(
                        ((Filter.FieldFilterElement) element).field,
                        getFilter(filterChain)
                );
        }
        return Filters.ALL;
    }

    @Override
    @SuppressWarnings("unchecked")
    public ObjectFilter getObjectFilter(Filter filter) {
        return getObjectFilter(new ArrayDeque<>(filter.filterChain));
    }

    @Override
    @SuppressWarnings("unchecked")
    public org.dizitart.no2.Filter getFilter(Filter filter) {
        return getFilter(new ArrayDeque<>(filter.filterChain));
    }
}
