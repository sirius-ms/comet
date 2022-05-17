package de.unijena.bioinf.ChemistryBase.ms.utils;

import de.unijena.bioinf.ChemistryBase.ms.Deviation;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Optional;

public class MassMap<T> {

    private final HashMap<Long, Entry<T>> map;
    private final double blowupFactor;

    public MassMap(double blowupFactor) {
        this.blowupFactor = blowupFactor;
        this.map = new HashMap<>();
    }

    public boolean put(double key, T value) {
        final long k = (int)(key*blowupFactor);
        return store(key, k, value);
    }

    private boolean store(double rawKey, long k, T value) {
        int size= map.size();
        Entry<T> e = map.computeIfAbsent(k, (x) -> new Entry<>(rawKey, value));
        if (e.value!=value) {
            while (e.successor!=null) {
                e = e.successor;
                if (e.value==value) return false;
            }
            e.successor = new Entry<>(rawKey, value);
            return true;
        } else return map.size()>size;
    }

    public List<T> retrieveAll(double key, Deviation deviation) {
        return retrieveAll(key, deviation.absoluteFor(key));
    }
    public List<T> retrieveAll(double key, double deviation) {
        final ArrayList<T> results = new ArrayList<>();
        final double left = key-deviation, right = key+deviation;
        final long lk = (int)(left*blowupFactor), rk = (int)(right*blowupFactor);
        for (long i=lk; i <= rk; ++i) {
            retrieveFrom(results, left, right, i);
        }
        return results;
    }
    public Optional<T> retrieveClosest(double key, Deviation deviation) {
        return retrieveClosest(key, deviation.absoluteFor(key));
    }
    public Optional<T> retrieveClosest(double key, double deviation) {
        Entry<T> best = null;
        double bestDistance = Double.POSITIVE_INFINITY;
        final double left = key-deviation, right = key+deviation;
        final long lk = (int)(left*blowupFactor), rk = (int)(right*blowupFactor);
        for (long i=lk; i <= rk; ++i) {
            best = retrieveFrom(key, left, right, i, best, bestDistance);
        }
        return best==null ? Optional.empty() : Optional.of(best.value);
    }

    private Entry<T> retrieveFrom(double key, double left, double right, long i, Entry<T> best, double bestDistance) {
        Entry<T> e = map.get(i);
        while (e!=null) {
            if (e.key >= left && e.key <= right) {
                double d = Math.abs(key-e.key);
                if (d < bestDistance) {
                    bestDistance = d;
                    best = e;
                }
            }
            e = e.successor;
        }
        return best;
    }

    private void retrieveFrom(ArrayList<T> results, double left, double right, long lk) {
        Entry<T> e = map.get(lk);
        while (e!=null) {
            if (e.key >= left && e.key <= right) {
                results.add(e.value);
            }
            e = e.successor;
        }
    }

    private static class Entry<T> {
        private final double key;
        private final T value;
        private Entry<T> successor;

        public Entry(double key, T entry) {
            this.key = key;
            this.value = entry;
            this.successor = null;
        }
    }
}
