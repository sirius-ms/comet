/*
 *
 *  This file is part of the SIRIUS library for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2013-2020 Kai Dührkop, Markus Fleischauer, Marcus Ludwig, Martin A. Hoffman, Fleming Kretschmer and Sebastian Böcker,
 *  Chair of Bioinformatics, Friedrich-Schiller University.
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

package de.unijena.bioinf.ms.middleware.controller;

import de.unijena.bioinf.ms.properties.PropertyManager;
import org.springframework.http.MediaType;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.RestController;

@RestController
public class VersionInfoController {

    @RequestMapping(value = "/version.json", method = RequestMethod.GET, produces = MediaType.APPLICATION_JSON_VALUE)
    public String getVersionInfo() {
        return "{" +
                "\"nightsky_api_version\": \""+PropertyManager.getProperty("de.unijena.bioinf.siriusNightsky.version")+"\"" +
                ", \"sirius_version\": \"" +  PropertyManager.getProperty("de.unijena.bioinf.siriusFrontend.version") + "\"" +
                ", \"sirius_lib_version\": \"" + PropertyManager.getProperty("de.unijena.bioinf.sirius.version") + "\"" +
                ", \"fingerid_version\": \"" + PropertyManager.getProperty("de.unijena.bioinf.fingerid.version") + "\"" +
                "}";
    }
}
