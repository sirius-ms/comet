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

import com.fasterxml.jackson.annotation.JsonInclude;
import de.unijena.bioinf.ChemistryBase.utils.DataQuality;
import de.unijena.bioinf.ms.persistence.model.core.QualityReport;
import io.swagger.v3.oas.annotations.media.Schema;
import lombok.Builder;
import lombok.Getter;
import lombok.Setter;

import java.util.LinkedHashMap;

@Getter
@Setter
@Builder
@JsonInclude(JsonInclude.Include.NON_NULL)
@Schema(name = "AlignedFeatureQualityExperimental",
        description = "EXPERIMENTAL: This schema is experimental and may be changed (or even removed) without notice until it is declared stable.")
public class AlignedFeatureQuality {
    /**
     * Id of the feature (aligned over runs) this quality information belongs to.
     */
    @Schema(nullable = false, requiredMode = Schema.RequiredMode.REQUIRED)
    protected String alignedFeatureId;

    /**
     * Overall Quality
     */
    @Schema(nullable = false, requiredMode = Schema.RequiredMode.REQUIRED)
    private DataQuality overallQuality;
    /**
     * Contains all pre-computation quality information that belong to
     * this feature (aligned over runs), such as information about the quality of the peak shape, MS2 spectrum etc.,
     */
    @Schema(nullable = false, requiredMode = Schema.RequiredMode.REQUIRED)
    private LinkedHashMap<String, QualityReport.Category> categories;
}
