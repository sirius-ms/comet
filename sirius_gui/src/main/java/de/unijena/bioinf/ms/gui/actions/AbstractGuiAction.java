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

package de.unijena.bioinf.ms.gui.actions;

import de.unijena.bioinf.ms.gui.SiriusGui;
import de.unijena.bioinf.ms.gui.mainframe.MainFrame;

import javax.swing.*;
import java.beans.PropertyChangeListener;

public abstract class AbstractGuiAction extends AbstractAction {
    protected final MainFrame mainFrame;
    protected final SiriusGui gui;

    public AbstractGuiAction(SiriusGui gui) {
        this.gui = gui;
        this.mainFrame = this.gui.getMainFrame();
    }

    public AbstractGuiAction(String name, SiriusGui gui) {
        super(name);
        this.gui = gui;
        this.mainFrame = this.gui.getMainFrame();
    }

    public AbstractGuiAction(String name, Icon icon, SiriusGui gui) {
        super(name, icon);
        this.gui = gui;
        this.mainFrame = this.gui.getMainFrame();
    }

    /**
     * For cleaning up resources such as listeners;
     */
    public void destroy(){
        if (this instanceof PropertyChangeListener c)
            this.gui.getConnectionMonitor().removePropertyChangeListener(c);
    }
}
