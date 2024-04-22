package de.unijena.bioinf.lcms.adducts;

import de.unijena.bioinf.ms.persistence.model.core.feature.AlignedFeatures;
import de.unijena.bioinf.ms.persistence.model.core.trace.MergedTrace;
import de.unijena.bioinf.ms.persistence.model.core.trace.SourceTrace;
import de.unijena.bioinf.ms.persistence.model.core.trace.TraceRef;
import it.unimi.dsi.fastutil.Pair;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.longs.Long2DoubleMap;
import it.unimi.dsi.fastutil.longs.Long2ObjectMap;
import it.unimi.dsi.fastutil.longs.Long2ObjectOpenHashMap;

import java.io.IOException;
import java.util.List;
import java.util.Optional;

public interface TraceProvider {

    public Optional<MergedTrace> getMergeTrace(AlignedFeatures feature) throws IOException;

    public Long2ObjectMap<SourceTrace> getSourceTraces(AlignedFeatures features);

    public Long2DoubleMap getIntensities(AlignedFeatures features);

    public Optional<Pair<TraceRef, SourceTrace>> getSourceTrace(AlignedFeatures features, long runId);

}
