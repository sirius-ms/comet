package de.unijena.bioinf.ChemistryBase.ms;

import de.unijena.bioinf.ms.annotations.Ms2ExperimentAnnotation;
import de.unijena.bioinf.ms.properties.ParameterConfig;

public abstract class ConfigAnnotation implements Ms2ExperimentAnnotation {
    // this are the the configs read from this file (@ParameterConfig)
    public final ParameterConfig config;

    protected ConfigAnnotation(ParameterConfig config) {
        this.config = config;
    }

}
