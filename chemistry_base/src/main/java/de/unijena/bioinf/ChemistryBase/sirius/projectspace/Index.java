package de.unijena.bioinf.ChemistryBase.sirius.projectspace;

import de.unijena.bioinf.ChemistryBase.ms.Ms2ExperimentAnnotation;

/**
 * Is used to order experiments
 */
public class Index implements Ms2ExperimentAnnotation {

    public static final Index NO_INDEX = new Index(-1);

    public final int index;

    public Index(int index) {
        this.index = index;
    }
}
