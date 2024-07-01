/*
 *
 *  This file is part of the SIRIUS library for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2013-2020 Kai Dührkop, Markus Fleischauer, Marcus Ludwig, Martin A. Hoffman, Fleming Kretschmer and Sebastian Böcker,
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

package de.unijena.bioinf.ms.middleware.model.features;

import de.unijena.bioinf.ms.middleware.model.spectra.BasicSpectrum;
import io.swagger.v3.oas.annotations.media.Schema;
import lombok.Builder;
import lombok.Getter;

import java.util.List;
import java.util.Set;

@Getter
@Builder
public class FeatureImport {
    @Schema(nullable = true, requiredMode = Schema.RequiredMode.NOT_REQUIRED)
    protected String name;

    /**
     * Externally provided FeatureId (by some preprocessing tool). This FeatureId is NOT used by SIRIUS but is stored to ease mapping information back to the source.
     */
    @Schema(nullable = true, requiredMode = Schema.RequiredMode.NOT_REQUIRED)
    protected String externalFeatureId;

    @Schema(nullable = false, requiredMode = Schema.RequiredMode.REQUIRED)
    protected Double ionMass;

    @Schema(nullable = false, requiredMode = Schema.RequiredMode.REQUIRED)
    protected int charge;

    /**
     * Detected adducts of this feature. Can be NULL or empty if no adducts are known.
     */
    @Schema(nullable = true, requiredMode = Schema.RequiredMode.NOT_REQUIRED)
    protected Set<String> detectedAdducts;

    @Schema(nullable = true)
    protected Double rtStartSeconds;
    @Schema(nullable = true)
    protected Double rtEndSeconds;

    /**
     * Mass Spec data of this feature (input data)
     */
    @Schema(nullable = true)
    protected BasicSpectrum mergedMs1;
    @Schema(requiredMode = Schema.RequiredMode.REQUIRED)
    protected List<BasicSpectrum> ms1Spectra;
    @Schema(requiredMode = Schema.RequiredMode.REQUIRED)
    protected List<BasicSpectrum> ms2Spectra;
}
