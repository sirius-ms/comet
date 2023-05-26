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

package de.unijena.bioinf.storage.db.nosql;

import java.io.Closeable;
import java.io.IOException;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

public interface Database<DocType> extends Closeable, AutoCloseable {

    enum SortOrder {
        ASCENDING, DESCENDING
    }

    //region CRUD operations

    <T> int insert(T object) throws IOException;

    <T> int insertAll(Iterable<T> objects) throws IOException;

    int insert(String collectionName, DocType document) throws IOException;

    int insertAll(String collectionName, Iterable<DocType> documents) throws IOException;

    <T> int upsert(T object) throws IOException;

    <T> int upsertAll(Iterable<T> objects) throws IOException;

    int upsert(String collectionName, DocType document) throws IOException;

    int upsertAll(String collectionName, Iterable<DocType> documents) throws IOException;

    <T> T getById(long id, Class<T> clazz) throws IOException;

    DocType getById(String collectionName, long id) throws IOException;

    <T> Iterable<T> find(Filter filter, Class<T> clazz) throws IOException;

    <T> Iterable<T> find(Filter filter, Class<T> clazz, int offset, int pageSize) throws IOException;

    <T> Iterable<T> find(Filter filter, Class<T> clazz, String sortField, SortOrder sortOrder) throws IOException;

    <T> Iterable<T> find(Filter filter, Class<T> clazz, int offset, int pageSize, String sortField, SortOrder sortOrder) throws IOException;

    Iterable<DocType> find(String collectionName, Filter filter) throws IOException;

    Iterable<DocType> find(String collectionName, Filter filter, int offset, int pageSize) throws IOException;

    Iterable<DocType> find(String collectionName, Filter filter, String sortField, SortOrder sortOrder) throws IOException;

    Iterable<DocType> find(String collectionName, Filter filter, int offset, int pageSize, String sortField, SortOrder sortOrder) throws IOException;

    <T> Iterable<T> findAll(Class<T> clazz) throws IOException;

    <T> Iterable<T> findAll(Class<T> clazz, int offset, int pageSize) throws IOException;

    <T> Iterable<T> findAll(Class<T> clazz, String sortField, SortOrder sortOrder) throws IOException;

    <T> Iterable<T> findAll(Class<T> clazz, int offset, int pageSize, String sortField, SortOrder sortOrder) throws IOException;

    Iterable<DocType> findAll(String collectionName) throws IOException;

    Iterable<DocType> findAll(String collectionName, int offset, int pageSize) throws IOException;

    Iterable<DocType> findAll(String collectionName, String sortField, SortOrder sortOrder) throws IOException;

    Iterable<DocType> findAll(String collectionName, int offset, int pageSize, String sortField, SortOrder sortOrder) throws IOException;

    <P, C> Iterable<P> joinAllChildren(Class<C> childClass, Iterable<P> parents, String localField, String foreignField, String targetField) throws IOException;

    <P, C> Iterable<P> joinChildren(Class<C> childClass, Filter childFilter, Iterable<P> parents, String localField, String foreignField, String targetField) throws IOException;

    <T, P, C> Iterable<T> joinAllChildren(Class<T> targetClass, Class<C> childClass, Iterable<P> parents, String localField, String foreignField, String targetField) throws IOException;

    <T, P, C> Iterable<T> joinChildren(Class<T> targetClass, Class<C> childClass, Filter childFilter, Iterable<P> parents, String localField, String foreignField, String targetField) throws IOException;

    Iterable<DocType> joinAllChildren(String childCollection, Iterable<DocType> parents, String localField, String foreignField, String targetField) throws IOException;

    Iterable<DocType> joinChildren(String childCollectionName, Filter childFilter, Iterable<DocType> parents, String localField, String foreignField, String targetField) throws IOException;

    <T> int count(Filter filter, Class<T> clazz) throws IOException;

    <T> int count(Filter filter, Class<T> clazz, int offset, int pageSize) throws IOException;

    <T> int countAll(Class<T> clazz) throws IOException;

    int count(String collectionName, Filter filter) throws IOException;

    int count(String collectionName, Filter filter, int offset, int pageSize) throws IOException;

    int countAll(String collectionName) throws IOException;

    <T> int remove(T object) throws IOException;

    <T> int removeAll(Iterable<T> objects) throws IOException;

    <T> int removeAll(Filter filter, Class<T> clazz) throws IOException;

    int remove(String collectionName, DocType document) throws IOException;

    int removeAll(String collectionName, Iterable<DocType> documents) throws IOException;

    int removeAll(String collectionName, Filter filter) throws IOException;

    //endregion
    <ObjectFilterType> ObjectFilterType getObjectFilter(Filter filter);

    <FilterType> FilterType getFilter(Filter filter);


    // region CRUD Streaming API
    default <T> Stream<T> findStr(Filter filter, Class<T> clazz) throws IOException {
        return StreamSupport.stream(find(filter, clazz).spliterator(), false);
    }

    default <T> Stream<T> findStr(Filter filter, Class<T> clazz, int offset, int pageSize) throws IOException {
        return StreamSupport.stream(find(filter, clazz, offset, pageSize).spliterator(), false);
    }

    default <T> Stream<T> findStr(Filter filter, Class<T> clazz, String sortField, SortOrder sortOrder) throws IOException {
        return StreamSupport.stream(find(filter, clazz, sortField, sortOrder).spliterator(), false);
    }

    default <T> Stream<T> findStr(Filter filter, Class<T> clazz, int offset, int pageSize, String sortField, SortOrder sortOrder) throws IOException {
        return StreamSupport.stream(find(filter, clazz, offset, pageSize, sortField, sortOrder).spliterator(), false);
    }

    default Stream<DocType> findStr(String collectionName, Filter filter) throws IOException {
        return StreamSupport.stream(find(collectionName, filter).spliterator(), false);
    }

    default Stream<DocType> findStr(String collectionName, Filter filter, int offset, int pageSize) throws IOException {
        return StreamSupport.stream(find(collectionName, filter, offset, pageSize).spliterator(), false);
    }

    default Stream<DocType> findStr(String collectionName, Filter filter, String sortField, SortOrder sortOrder) throws IOException {
        return StreamSupport.stream(find(collectionName, filter, sortField, sortOrder).spliterator(), false);
    }

    default Stream<DocType> findStr(String collectionName, Filter filter, int offset, int pageSize, String sortField, SortOrder sortOrder) throws IOException {
        return StreamSupport.stream(find(collectionName, filter, offset, pageSize, sortField, sortOrder).spliterator(), false);
    }

    default <T> Stream<T> findAllStr(Class<T> clazz) throws IOException {
        return StreamSupport.stream(findAll(clazz).spliterator(), false);
    }

    default <T> Stream<T> findAllStr(Class<T> clazz, int offset, int pageSize) throws IOException {
        return StreamSupport.stream(findAll(clazz, offset, pageSize).spliterator(), false);
    }

    default <T> Stream<T> findAllStr(Class<T> clazz, String sortField, SortOrder sortOrder) throws IOException {
        return StreamSupport.stream(findAll(clazz, sortField, sortOrder).spliterator(), false);
    }

    default <T> Stream<T> findAllStr(Class<T> clazz, int offset, int pageSize, String sortField, SortOrder sortOrder) throws IOException {
        return StreamSupport.stream(findAll(clazz, offset, pageSize, sortField, sortOrder).spliterator(), false);
    }

    default Stream<DocType> findAllStr(String collectionName) throws IOException {
        return StreamSupport.stream(findAll(collectionName).spliterator(), false);
    }

    default Stream<DocType> findAllStr(String collectionName, int offset, int pageSize) throws IOException {
        return StreamSupport.stream(findAll(collectionName, offset, pageSize).spliterator(), false);
    }

    default Stream<DocType> findAllStr(String collectionName, String sortField, SortOrder sortOrder) throws IOException {
        return StreamSupport.stream(findAll(collectionName, sortField, sortOrder).spliterator(), false);
    }

    default Stream<DocType> findAllStr(String collectionName, int offset, int pageSize, String sortField, SortOrder sortOrder) throws IOException {
        return StreamSupport.stream(findAll(collectionName, offset, pageSize, sortField, sortOrder).spliterator(), false);
    }

    default <T, P, C> Stream<T> joinAllChildrenStr(Class<T> targetClass, Class<C> childClass, Iterable<P> parents, String localField, String foreignField, String targetField) throws IOException {
        return StreamSupport.stream(joinAllChildren(targetClass, childClass, parents, localField, foreignField, targetField).spliterator(), false);
    }

    default <T, P, C> Stream<T> joinChildrenStr(Class<T> clazz, Class<C> childClass, Filter childFilter, Iterable<P> parents, String localField, String foreignField, String targetField) throws IOException {
        return StreamSupport.stream(joinChildren(clazz, childClass, childFilter, parents, localField, foreignField, targetField).spliterator(), false);
    }

    default Stream<DocType> joinAllChildrenStr(String childCollectionName, Iterable<DocType> parents, String localField, String foreignField, String targetField) throws IOException  {
        return StreamSupport.stream(joinAllChildren(childCollectionName, parents, localField, foreignField, targetField).spliterator(), false);
    }

    default Stream<DocType> joinChildrenStr(String childCollectionName, Filter childFilter, Iterable<DocType> parents, String localField, String foreignField, String targetField) throws IOException {
        return StreamSupport.stream(joinChildren(childCollectionName, childFilter, parents, localField, foreignField, targetField).spliterator(), false);
    }
    //endregion
}
