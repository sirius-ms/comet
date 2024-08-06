package de.unijena.bioinf.fragmenter;

import de.unijena.bioinf.ms.annotations.Ms2ExperimentAnnotation;
import de.unijena.bioinf.ms.properties.DefaultInstanceProvider;
import de.unijena.bioinf.ms.properties.DefaultProperty;

public class RankWithEpimetheus implements Ms2ExperimentAnnotation {

    public static final RankWithEpimetheus TRUE = new RankWithEpimetheus(true);
    public static final RankWithEpimetheus FALSE = new RankWithEpimetheus(false);

    public final boolean value;

    private RankWithEpimetheus(boolean value){
        this.value = value;
    }
    @DefaultInstanceProvider
    public static RankWithEpimetheus newInstance(@DefaultProperty boolean value){
        return value ? TRUE : FALSE;
    }
}
