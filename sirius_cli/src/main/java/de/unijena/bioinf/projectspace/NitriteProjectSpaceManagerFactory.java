/*
 *
 *  This file is part of the SIRIUS library for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2023 Bright Giant GmbH
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
 *  You should have received a copy of the GNU General Public License along with SIRIUS.
 *  If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>
 */

package de.unijena.bioinf.projectspace;

import com.github.f4b6a3.tsid.TsidCreator;
import de.unijena.bioinf.ChemistryBase.utils.FileUtils;
import de.unijena.bioinf.ms.persistence.storage.nitrite.NitriteSirirusProject;
import lombok.extern.slf4j.Slf4j;
import org.jetbrains.annotations.Nullable;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import static de.unijena.bioinf.ms.persistence.storage.SiriusProjectDocumentDatabase.SIRIUS_PROJECT_SUFFIX;

@Slf4j
public class NitriteProjectSpaceManagerFactory implements ProjectSpaceManagerFactory<NoSQLProjectSpaceManager> {

    @Override
    public NoSQLProjectSpaceManager createOrOpen(@Nullable Path projectLocation) throws IOException {

        if (projectLocation == null) {
            projectLocation = FileUtils.newTempFile("sirius-tmp-project-" + TsidCreator.getTsid(), SIRIUS_PROJECT_SUFFIX);
            log.warn("No unique output location found. Writing output to Temporary folder: " + projectLocation.toString());
            if (Files.exists(projectLocation)) {
                throw new IOException("Could not create new Project '" + projectLocation + "' because it already exists");
            }
        }

        return new NoSQLProjectSpaceManager(new NitriteSirirusProject(projectLocation));
    }

}
