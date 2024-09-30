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

import de.unijena.bioinf.ms.middleware.model.tags.TagCategory;
import de.unijena.bioinf.ms.middleware.service.projects.ProjectsProvider;
import io.swagger.v3.oas.annotations.tags.Tag;
import jakarta.validation.Valid;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.MediaType;
import org.springframework.web.bind.annotation.*;

import java.util.List;

@RestController
@RequestMapping(value = "/api/projects/{projectId}/categories")
@Tag(name = "Tag categories", description = "This API allows accessing tag categories.")
public class CategoryController {

    final protected ProjectsProvider<?> projectsProvider;

    @Autowired
    public CategoryController(ProjectsProvider<?> projectsProvider) {
        this.projectsProvider = projectsProvider;
    }

    /**
     * Get all tag categories in the given project-space.
     *
     * @param projectId project-space to read from.
     * @return Tag categories.
     */
    @GetMapping(produces = MediaType.APPLICATION_JSON_VALUE)
    public List<TagCategory> getCategories(@PathVariable String projectId) {
        return projectsProvider.getProjectOrThrow(projectId).findCategories();
    }

    /**
     * Get tag categories by type in the given project-space.
     *
     * @param projectId    project-space to read from.
     * @param categoryType name of the category
     * @return Tag categories.
     */
    @GetMapping(value = "/type/{categoryType}", produces = MediaType.APPLICATION_JSON_VALUE)
    public List<TagCategory> getCategoriesByType(@PathVariable String projectId, @PathVariable String categoryType) {
        return projectsProvider.getProjectOrThrow(projectId).findCategoriesByType(categoryType);
    }

    /**
     * Get tag category by name in the given project-space.
     *
     * @param projectId    project-space to read from.
     * @param categoryName name of the category
     * @return Tag categories.
     */
    @GetMapping(value = "/name/{categoryName}", produces = MediaType.APPLICATION_JSON_VALUE)
    public TagCategory getCategoryByName(@PathVariable String projectId, @PathVariable String categoryName) {
        return projectsProvider.getProjectOrThrow(projectId).findCategoryByName(categoryName);
    }

    /**
     * Add tag category to the project. Category name must not exist in the project.
     *
     * @param projectId  project-space to add to.
     * @param categories the tag categories to be added
     * @return the tag categories that have been added
     */
    @PostMapping(value = "/add", produces = MediaType.APPLICATION_JSON_VALUE)
    public List<TagCategory> addCategories(@PathVariable String projectId, @Valid @RequestBody List<TagCategory> categories) {
        return projectsProvider.getProjectOrThrow(projectId).addCategories(categories);
    }

    /**
     * Delete tag categories with the given names from the specified project-space.
     *
     * @param projectId     project-space to delete from.
     * @param categoryNames names of categories to delete.
     */
    @PutMapping(value = "/delete")
    public void deleteCategories(@PathVariable String projectId, @RequestBody List<String> categoryNames) {
        projectsProvider.getProjectOrThrow(projectId).deleteCategories(categoryNames);
    }

}
