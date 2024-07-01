package de.unijena.bioinf.ms.gui.utils;/*
 *
 *  This file is part of the SIRIUS library for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2013-2021 Kai Dührkop, Markus Fleischauer, Marcus Ludwig, Martin A. Hoffman and Sebastian Böcker,
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
 *  You should have received a copy of the GNU General Public License along with SIRIUS. If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>
 */

import ca.odell.glazedlists.matchers.AbstractMatcherEditor;
import de.unijena.bioinf.projectspace.InstanceBean;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;

public class CompoundFilterMatcherEditor extends AbstractMatcherEditor<InstanceBean> implements PropertyChangeListener {
    final CompoundFilterMatcher matchter;

    public CompoundFilterMatcherEditor(CompoundFilterMatcher matchter) {
        this.matchter = matchter;
        this.matchter.filterModel.addPropertyChangeListener("filterUpdateCompleted", this);
    }

    @Override
    public void propertyChange(PropertyChangeEvent evt) {
        if (evt.getSource() == matchter.filterModel)
            fireChanged(matchter);
    }
}
