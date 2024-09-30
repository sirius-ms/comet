/*
 *
 *  This file is part of the SIRIUS library for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2024 Bright Giant GmbH
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

package de.unijena.bioinf.ms.middleware.controller;

import de.unijena.bioinf.ms.middleware.model.tags.TagFilter;
import de.unijena.bioinf.ms.middleware.service.projects.ProjectsProvider;
import io.swagger.v3.oas.annotations.tags.Tag;
import jakarta.validation.Valid;
import org.jetbrains.annotations.Nullable;
import org.springdoc.core.annotations.ParameterObject;
import org.springframework.data.domain.Page;
import org.springframework.data.domain.Pageable;
import org.springframework.http.MediaType;
import org.springframework.web.bind.annotation.*;

import java.util.EnumSet;
import java.util.List;

public abstract class TagController<T, O extends Enum<O>> {

    final protected ProjectsProvider<?> projectsProvider;

    protected TagController(ProjectsProvider<?> projectsProvider) {
        this.projectsProvider = projectsProvider;
    }

    protected abstract Class<T> getTaggable();

    @PostMapping(value = "/tags/tagged/{categoryName}", produces = MediaType.APPLICATION_JSON_VALUE)
    public Page<T> objectsByTag(@PathVariable String projectId,
                                @PathVariable String categoryName,
                                @Nullable @RequestBody TagFilter filter,
                                @ParameterObject Pageable pageable,
                                @RequestParam(defaultValue = "") EnumSet<O> optFields
    ) {
        return projectsProvider.getProjectOrThrow(projectId).findObjectsByTag(getTaggable(), categoryName, filter, pageable, optFields);
    }

    /**
     * Add tags to a run in the project. Tags with the same category name will be overwritten.
     *
     * @param projectId  project-space to add to.
     * @param objectId   object to tag.
     * @param tags       tags to add.
     * @return the tags that have been added
     */
    @PostMapping(value = "/tags/add/{objectId}", produces = MediaType.APPLICATION_JSON_VALUE)
    public List<Tag> addTags(@PathVariable String projectId, @PathVariable String objectId, @Valid @RequestBody List<Tag> tags) {
        return projectsProvider.getProjectOrThrow(projectId).addTagsToObject(getTaggable(), objectId, tags);
    }

    /**
     * Delete tags with the given IDs from the specified project-space.
     *
     * @param projectId     project-space to delete from.
     * @param objectId      object to delete tags from.
     * @param categoryNames Category names of the tags to delete.
     */
    @PutMapping(value = "/tags/delete/{objectId}")
    public void deleteTags(@PathVariable String projectId, @PathVariable String objectId, @RequestBody List<String> categoryNames) {
        projectsProvider.getProjectOrThrow(projectId).deleteTagsFromObject(getTaggable(), objectId, categoryNames);
    }

}
