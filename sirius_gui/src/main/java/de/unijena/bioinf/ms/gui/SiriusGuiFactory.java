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


package de.unijena.bioinf.ms.gui;

import de.unijena.bioinf.ms.gui.net.ConnectionMonitor;
import de.unijena.bioinf.ms.nightsky.sdk.NightSkyClient;
import org.jetbrains.annotations.NotNull;

/**
 * This is usually a singleton instance class that holds all active gui instances and there shared infrastructure
 */
public final class SiriusGuiFactory {
    private volatile NightSkyClient nightSkyClient;
    private volatile ConnectionMonitor connectionMonitor;


    public SiriusGuiFactory() {
        this(null, null);
    }
    public SiriusGuiFactory(NightSkyClient nightSkyClient, ConnectionMonitor connectionMonitor) {
        this.nightSkyClient = nightSkyClient;
        this.connectionMonitor = connectionMonitor;
    }

    public SiriusGui newGui(@NotNull String projectId){
        init();
        return new SiriusGui(projectId, nightSkyClient, connectionMonitor);
    }

    private void init(){
        if (nightSkyClient == null){
            synchronized (this){
                if (nightSkyClient == null)
                    nightSkyClient = new NightSkyClient();
            }
        }

        if (connectionMonitor == null){
            synchronized (this){
                if (connectionMonitor == null)
                    connectionMonitor = new ConnectionMonitor(nightSkyClient);
            }
        }
    }

    public void shutdowm(){
        connectionMonitor.close();
    }
}
