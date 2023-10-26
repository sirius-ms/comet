/*
 *
 *  This file is part of the SIRIUS library for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2013-2020 Kai Dührkop, Markus Fleischauer, Marcus Ludwig, Martin A. Hoffman, Fleming Kretschmer and Sebastian Böcker,
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
 *  You should have received a copy of the GNU Lesser General Public License along with SIRIUS. If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>
 */

package de.unijena.bioinf.ms.middleware.model.features.annotations;

import de.unijena.bioinf.ChemistryBase.fp.NPCFingerprintVersion;

public class CanopusLevels {
    public static final String[] classyFireLevelNames = new String[]{
            "Kingdom", "Superclass", "Class", "Subclass"
    };

    public static String getNPCLevelName(int level){
        if (level > 2)
            throw new IllegalArgumentException("NPC level must be between 0 and 2");
        return NPCFingerprintVersion.NPCLevel.values()[level].name();

    }
    public static String getClassyFireLevelName(int level){
        if (level > classyFireLevelNames.length)
            return "Level-" + level;
        return classyFireLevelNames[level-1];
    }
}
