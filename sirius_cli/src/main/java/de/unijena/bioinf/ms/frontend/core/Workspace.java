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

package de.unijena.bioinf.ms.frontend.core;

import de.unijena.bioinf.ms.properties.PropertyManager;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;

public class Workspace {
    public static final Path WORKSPACE;

    public static final Path loggingPropFile;
    public static final Path siriusPropsFile;
    public static final Path customProfileFile;
    public static final Path versionFile;

    public static final Path runConfigDir;

    static {
        try {
            System.setProperty("de.unijena.bioinf.ms.propertyLocations", "sirius_frontend.build.properties");

            final String version = PropertyManager.getProperty("de.unijena.bioinf.siriusFrontend.version");
            final String[] versionParts = version.split("[.]");
            //#################### start init workspace ################################
            final String home = System.getProperty("user.home");
            String defaultFolderName = PropertyManager.getProperty("de.unijena.bioinf.sirius.ws.default.name", null, ".sirius");
            if (versionParts != null && versionParts.length > 1)
                defaultFolderName = defaultFolderName + "-" + versionParts[0] + "." + versionParts[1];

            final Path DEFAULT_WORKSPACE = Paths.get(home).resolve(defaultFolderName);
            final Map<String, String> env = System.getenv();

            String ws = System.getProperty("de.unijena.bioinf.sirius.ws.location"); //todo dirty hack make cli parsing work correclty
            if (ws == null)
                ws = PropertyManager.getProperty("de.unijena.bioinf.sirius.ws.location");
            if (ws == null)
                ws = env.get("SIRIUS_WORKSPACE");
            if (ws != null) {
                Path wsDir = Paths.get(ws);
                if (Files.isDirectory(wsDir)) {
                    WORKSPACE = wsDir;
                } else if (Files.notExists(wsDir)) {
                    try {
                        Files.createDirectories(wsDir);
                    } catch (IOException e) {
                        System.err.println("Could not create Workspace set in environment variable! Falling back to default Workspace - " + DEFAULT_WORKSPACE.toString());
                        e.printStackTrace();
                        wsDir = DEFAULT_WORKSPACE;
                    } finally {
                        WORKSPACE = wsDir;
                    }
                } else {
                    System.err.println("WARNING: " + wsDir.toString() + " is not a directory! Falling back to default Workspace - " + DEFAULT_WORKSPACE);
                    WORKSPACE = DEFAULT_WORKSPACE;
                }
            } else {
                WORKSPACE = DEFAULT_WORKSPACE;
            }

            if (Files.notExists(WORKSPACE)) {
                try {
                    Files.createDirectories(WORKSPACE);
                } catch (IOException e) {
                    System.err.println("Could NOT create Workspace");
                    e.printStackTrace();
                    System.exit(1);
                }
            }

            loggingPropFile = WORKSPACE.resolve("logging.properties");
            siriusPropsFile = WORKSPACE.resolve("sirius.properties");
            customProfileFile = WORKSPACE.resolve("custom.config");
            versionFile = WORKSPACE.resolve("version");
            runConfigDir = WORKSPACE.resolve("run-configs");

            Files.createDirectories(runConfigDir);
        } catch (Exception e) {
            System.err.println("Workspace Core STATIC Block Error!");
            e.printStackTrace(System.err);
            throw new RuntimeException(e);
        }
    }
}
