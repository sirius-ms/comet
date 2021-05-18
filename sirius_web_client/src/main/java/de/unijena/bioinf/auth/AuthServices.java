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

package de.unijena.bioinf.auth;

import de.unijena.bioinf.auth.auth0.Auth0Api;
import de.unijena.bioinf.babelms.utils.Base64;
import de.unijena.bioinf.ms.properties.PropertyManager;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.net.MalformedURLException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

public class AuthServices {
    private AuthServices() {
    }

    public static AuthService createDefault(Path refreshTokenFile) throws MalformedURLException {
        String rToken = null;
        try {
            if (Files.isExecutable(refreshTokenFile))
                rToken = readRefreshToken(refreshTokenFile);
        } catch (IOException e) {
            LoggerFactory.getLogger(AuthServices.class).warn("Could not read refresh token from file! (re)login might be needed!", e);
        }

        Auth0Api api = new Auth0Api(PropertyManager.getProperty("de.unijena.bioinf.sirius.security.authServer"));
        return new AuthService(rToken, api);
    }

    public static void writeRefreshToken(String refreshToken, Path refreshTokenFile) throws IOException {
        Files.write(refreshTokenFile, Base64.encodeBytesToBytes(refreshToken.getBytes(StandardCharsets.UTF_8)));
    }

    public static String readRefreshToken(Path refreshTokenFile) throws IOException {
        return new String(Base64.decode(Files.readAllBytes(refreshTokenFile)), StandardCharsets.UTF_8);
    }
}
