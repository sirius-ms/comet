package de.unijena.bioinf.lcms.adducts;

import de.unijena.bioinf.ms.persistence.model.core.feature.AlignedFeatures;
import de.unijena.bioinf.ms.persistence.model.core.feature.Feature;
import de.unijena.bioinf.ms.persistence.model.core.trace.MergedTrace;
import de.unijena.bioinf.ms.persistence.model.core.trace.SourceTrace;
import de.unijena.bioinf.ms.persistence.model.core.trace.TraceRef;
import de.unijena.bioinf.ms.persistence.storage.MsProjectDocumentDatabase;
import de.unijena.bioinf.storage.db.nosql.Filter;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.longs.Long2DoubleMap;
import it.unimi.dsi.fastutil.longs.Long2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.longs.Long2ObjectMap;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;

import java.io.IOException;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;

public class ProjectSpaceTraceProvider implements TraceProvider {

    protected final MsProjectDocumentDatabase<?> storage;

    public ProjectSpaceTraceProvider(MsProjectDocumentDatabase<?> storage) {
        this.storage = storage;
    }

    @Override
    public Optional<MergedTrace> getMergeTrace(AlignedFeatures feature) {
        return feature.getTraceRef().map(id-> {
            Iterator<MergedTrace> iter = null;
            try {
                iter = storage.getStorage().find(Filter.where("mergedTraceId").eq(id), MergedTrace.class).iterator();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
            if (iter.hasNext()) {
                return iter.next();
            } else return null;
        });
    }

    @Override
    public Long2ObjectMap<SourceTrace> getSourceTraces(AlignedFeatures features) {
        Long2ObjectOpenHashMap<SourceTrace> traces = new Long2ObjectOpenHashMap<>();
        for (Feature f : getFeatures(features)) {
            if (f.getTraceRef().isPresent()) {
                TraceRef traceRef = f.getTraceRef().get();
                Iterator<SourceTrace> iter = null;
                try {
                    iter = storage.getStorage().find(Filter.where("sourceTraceId").eq(traceRef.getTraceId()), SourceTrace.class).iterator();
                    if (iter.hasNext()) {
                        traces.put(f.getRunId(), iter.next());
                    }
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
        }
        return traces;
    }

    @Override
    public Long2DoubleMap getIntensities(AlignedFeatures features) {
        Long2DoubleMap map = new Long2DoubleOpenHashMap();
        for (Feature f : getFeatures(features)) {
            map.put(f.getRunId(), f.getApexIntensity());
        }
        return map;
    }

    private List<Feature> getFeatures(AlignedFeatures feature) {
        if (feature.getFeatures().isEmpty()) return Collections.emptyList();
        else return feature.getFeatures().get();
    }
}
