/*
 *  This file is part of the SIRIUS Software for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2013-2020 Kai Dührkop, Markus Fleischauer, Marcus Ludwig, Martin A. Hoffman, Fleming Kretschmer, Marvin Meusel and Sebastian Böcker,
 *  Chair of Bioinformatics, Friedrich-Schiller University.
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Affero General Public License
 *  as published by the Free Software Foundation; either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License along with SIRIUS.  If not, see <https://www.gnu.org/licenses/agpl-3.0.txt>
 */

package de.unijena.bioinf.ms.gui.molecular_formular;

import de.unijena.bioinf.ms.gui.table.SiriusTableFormat;
import de.unijena.bioinf.ms.nightsky.sdk.model.Deviation;
import de.unijena.bioinf.ms.nightsky.sdk.model.FormulaCandidate;
import de.unijena.bioinf.ms.nightsky.sdk.model.FormulaCandidateOptField;
import de.unijena.bioinf.ms.nightsky.sdk.model.LipidAnnotation;
import de.unijena.bioinf.projectspace.FormulaResultBean;
import org.jetbrains.annotations.NotNull;
import org.slf4j.LoggerFactory;

import java.util.OptionalInt;
import java.util.function.Function;

/**
 * Display issues in a tabular form.
 *
 * @author Markus Fleischauer (markus.fleischauer@gmail.com)
 */
public class SiriusResultTableFormat extends SiriusTableFormat<FormulaResultBean> {
    private static final int COL_COUNT = 11;

    protected SiriusResultTableFormat(Function<FormulaResultBean,Boolean> isBest) {
        super(isBest);
    }

    @Override
    public int highlightColumnIndex() {
        return COL_COUNT;
    }

    @Override
    public int getColumnCount() {
        return COL_COUNT;
    }

    @Override
    public String getColumnName(int column) {
        return switch (column) {
            case 0 -> "Molecular Formula";
            case 1 -> "Adduct";
            case 2 -> "Zodiac Score";
            case 3 -> "Sirius Score";
            case 4 -> "Isotope Score";
            case 5 -> "Tree Score";
            case 6 -> "Explained Peaks";
            case 7 -> "Total Explained Intensity";
            case 8 -> "Median Mass Error (ppm)";
            case 9 -> "Median Absolute Mass Error (ppm)";
            case 10 -> "Lipid Class";
            case 11 -> "Best";
            default -> throw new IllegalStateException();
        };
    }

    @Override
    public Object getColumnValue(FormulaResultBean result, int column) {
        return switch (column) {
            case 0 -> result.getMolecularFormula();
            case 1 -> result.getAdduct();
            case 2 -> result.getZodiacScore().orElse(Double.NaN);
            case 3 -> result.getSiriusScore().orElse(Double.NaN);
            case 4 -> result.getIsotopeScore().orElse(Double.NaN);
            case 5 -> result.getTreeScore().orElse(Double.NaN);
            case 6 -> result.getNumOfExplainedPeaks().stream().mapToDouble(v -> (double) v).findFirst().orElse(Double.NaN);
            case 7 -> result.getTotalExplainedIntensity().orElse(Double.NaN);
            case 8 -> result.getMedianMassDeviation().map(Deviation::getPpm).orElse(Double.NaN);
            case 9 -> result.getMedianMassDeviation().map(Deviation::getAbsolute).orElse(Double.NaN);
            case 10 -> result.getLipidAnnotation().map(LipidAnnotation::getLipidSpecies).orElse("N/A"); //N/A or better empty?
            case 11 -> isBest.apply(result);
            default -> throw new IllegalStateException();
        };
    }
}

