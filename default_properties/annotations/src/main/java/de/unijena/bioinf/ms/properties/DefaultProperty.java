/*
 *
 *  This file is part of the SIRIUS library for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2013-2020 Kai Dührkop, Markus Fleischauer, Marcus Ludwig, Martin A. Hoffman and Sebastian Böcker,
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
 *  You should have received a copy of the GNU General Public License along with SIRIUS. If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>
 */

package de.unijena.bioinf.ms.properties;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

@Retention(RetentionPolicy.RUNTIME)
@Target({ElementType.TYPE, ElementType.FIELD, ElementType.PARAMETER})
public @interface DefaultProperty {
    /**
     * Use this if you want the id for the given field or parameter to be
     * different from the field name. Note that the field name will be
     * ignored if a Class has only a single field/parameter. In such
     * cases only the propertyParent value will be used.
     * */
    String propertyKey() default "";

    /**
     * Use this if you wont the ID/Name of the Property to be different from the corresponding Class name.
     * One Class/Parent may hold multiple fields.
     * Note that all propertyParent values in one Class need to be the same and that the
     * propertyParent name is used as a unique ID to identify the corresponding Class at runtime.
     * */
    String propertyParent() default "";
}
