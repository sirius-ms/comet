/*
 *
 *  This file is part of the SIRIUS library for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2013-2020 Kai Dührkop, Markus Fleischauer, Marcus Ludwig, Martin A. Hoffman and Sebastian Böcker,
 *  Chair of Bioinformatics, Friedrich-Schilller University.
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
 *  You should have received a copy of the GNU General Public License along with SIRIUS. If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>
 */

package de.unijena.bioinf.ms.rest.model.worker;


import com.fasterxml.jackson.annotation.JsonAutoDetect;
import io.swagger.v3.oas.annotations.media.Schema;
import org.jetbrains.annotations.NotNull;

import java.time.Instant;
import java.util.*;
import java.util.function.Consumer;
import java.util.function.Predicate;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static io.swagger.v3.oas.annotations.media.Schema.RequiredMode.REQUIRED;

@JsonAutoDetect(fieldVisibility = JsonAutoDetect.Visibility.ANY, getterVisibility = JsonAutoDetect.Visibility.NONE, setterVisibility = JsonAutoDetect.Visibility.NONE, isGetterVisibility = JsonAutoDetect.Visibility.NONE)
public class WorkerList {
    @Schema(nullable = false, requiredMode = REQUIRED)
    private int pendingJobs = Integer.MIN_VALUE;
    @Schema(nullable = false, requiredMode = REQUIRED)
    private final ArrayList<WorkerInfo> workerList;

    public WorkerList(int initialCapacity) {
        workerList = new ArrayList<>(initialCapacity);
    }

    public WorkerList(Collection<? extends WorkerInfo> c, int pendingJobs) {
        this(c);
        this.pendingJobs = pendingJobs;
    }

    public WorkerList() {
        workerList = new ArrayList<>();
    }


    public WorkerList(Collection<? extends WorkerInfo> c) {
        workerList = new ArrayList<>(c);
    }

    public WorkerList(int initialCapacity, int pendingJobs) {
        this(initialCapacity);
        this.pendingJobs = pendingJobs;
    }

    public int getPendingJobs() {
        return pendingJobs;
    }

    public void setPendingJobs(int pendingJobs) {
        this.pendingJobs = pendingJobs;
    }

    public Stream<WorkerInfo> getWorkerActiveWithinAsStrean(Instant slot) {
        return workerList.stream().filter((w) -> w.isAlive(slot.toEpochMilli()));
    }

    public List<WorkerInfo> getWorkerActiveWithin(Instant slot) {
        return getWorkerActiveWithinAsStrean(slot).collect(Collectors.toList());
    }

    public long getNumWorkerActiveWithin(Instant slot) {
        return getWorkerActiveWithinAsStrean(slot).count();
    }

    public Set<WorkerWithCharge> getSupportedTypes() {
        return workerList.stream().map(WorkerInfo::asWorkerWithCharge).collect(Collectors.toSet());
    }

    public Set<WorkerWithCharge> getActiveSupportedTypes() {
        return getActiveSupportedTypes(Instant.ofEpochSecond(600/*10 min*/));
    }

    public Set<WorkerWithCharge> getActiveSupportedTypes(Instant slot) {
        return getWorkerActiveWithinAsStrean(slot).map(WorkerInfo::asWorkerWithCharge).collect(Collectors.toSet());
    }

    public boolean supportsAllPredictorTypes(Set<WorkerWithCharge> neededTypes, Instant activeWithin) {
        return getActiveSupportedTypes(activeWithin).containsAll(neededTypes);
    }

    public boolean supportsAllPredictorTypes(Set<WorkerWithCharge> neededTypes) {
        return getActiveSupportedTypes(Instant.ofEpochSecond(600/*10 min*/)).containsAll(neededTypes);
    }


    //region ArrayList Delegation
    public int size() {
        return workerList.size();
    }

    public boolean isEmpty() {
        return workerList.isEmpty();
    }

    public boolean contains(Object o) {
        return workerList.contains(o);
    }

    public int indexOf(Object o) {
        return workerList.indexOf(o);
    }

    public int lastIndexOf(Object o) {
        return workerList.lastIndexOf(o);
    }

    public WorkerInfo get(int index) {
        return workerList.get(index);
    }

    public WorkerInfo set(int index, WorkerInfo element) {
        return workerList.set(index, element);
    }

    public boolean add(WorkerInfo workerInfo) {
        return workerList.add(workerInfo);
    }

    public void add(int index, WorkerInfo element) {
        workerList.add(index, element);
    }

    public WorkerInfo remove(int index) {
        return workerList.remove(index);
    }

    public boolean remove(Object o) {
        return workerList.remove(o);
    }

    public void clear() {
        workerList.clear();
    }

    public boolean addAll(Collection<? extends WorkerInfo> c) {
        return workerList.addAll(c);
    }

    public boolean addAll(int index, Collection<? extends WorkerInfo> c) {
        return workerList.addAll(index, c);
    }

    public boolean removeAll(Collection<?> c) {
        return workerList.removeAll(c);
    }

    public boolean retainAll(Collection<?> c) {
        return workerList.retainAll(c);
    }

    @NotNull
    public ListIterator<WorkerInfo> listIterator(int index) {
        return workerList.listIterator(index);
    }

    @NotNull
    public ListIterator<WorkerInfo> listIterator() {
        return workerList.listIterator();
    }

    @NotNull
    public Iterator<WorkerInfo> iterator() {
        return workerList.iterator();
    }

    @NotNull
    public List<WorkerInfo> subList(int fromIndex, int toIndex) {
        return workerList.subList(fromIndex, toIndex);
    }

    public void forEach(Consumer<? super WorkerInfo> action) {
        workerList.forEach(action);
    }

    public Spliterator<WorkerInfo> spliterator() {
        return workerList.spliterator();
    }

    public boolean removeIf(Predicate<? super WorkerInfo> filter) {
        return workerList.removeIf(filter);
    }

    public void replaceAll(UnaryOperator<WorkerInfo> operator) {
        workerList.replaceAll(operator);
    }

    public boolean containsAll(Collection<?> c) {
        return workerList.containsAll(c);
    }

    public Stream<WorkerInfo> stream() {
        return workerList.stream();
    }

    public Stream<WorkerInfo> parallelStream() {
        return workerList.parallelStream();
    }
    //endregion
}
