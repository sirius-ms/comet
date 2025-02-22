/*
 *  This file is part of the SIRIUS libraries for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2024 Bright Giant GmbH
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with SIRIUS.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  NOTE: This class is auto generated by OpenAPI Generator (https://openapi-generator.tech).
 *  https://openapi-generator.tech
 *  Do not edit the class manually.
 */


package io.sirius.ms.sdk.model;

import java.util.Objects;
import java.util.Arrays;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonTypeName;
import com.fasterxml.jackson.annotation.JsonValue;
import io.sirius.ms.sdk.model.CompoundClasses;
import io.sirius.ms.sdk.model.ConfidenceMode;
import io.sirius.ms.sdk.model.FormulaCandidate;
import io.sirius.ms.sdk.model.StructureCandidateScored;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;
import com.fasterxml.jackson.annotation.JsonTypeName;

/**
 * Summary of the results of a feature (aligned over runs). Can be added to a AlignedFeature.  The different annotation fields within this summary object are null if the corresponding  feature does not contain the represented results. If fields are non-null  the corresponding result has been computed but might still be empty.
 */
@JsonPropertyOrder({
  FeatureAnnotations.JSON_PROPERTY_FORMULA_ANNOTATION,
  FeatureAnnotations.JSON_PROPERTY_STRUCTURE_ANNOTATION,
  FeatureAnnotations.JSON_PROPERTY_COMPOUND_CLASS_ANNOTATION,
  FeatureAnnotations.JSON_PROPERTY_CONFIDENCE_EXACT_MATCH,
  FeatureAnnotations.JSON_PROPERTY_CONFIDENCE_APPROX_MATCH,
  FeatureAnnotations.JSON_PROPERTY_EXPANSIVE_SEARCH_STATE,
  FeatureAnnotations.JSON_PROPERTY_SPECIFIED_DATABASES,
  FeatureAnnotations.JSON_PROPERTY_EXPANDED_DATABASES
})
@jakarta.annotation.Generated(value = "org.openapitools.codegen.languages.JavaClientCodegen", comments = "Generator version: 7.6.0")
public class FeatureAnnotations {
  public static final String JSON_PROPERTY_FORMULA_ANNOTATION = "formulaAnnotation";
  private FormulaCandidate formulaAnnotation;

  public static final String JSON_PROPERTY_STRUCTURE_ANNOTATION = "structureAnnotation";
  private StructureCandidateScored structureAnnotation;

  public static final String JSON_PROPERTY_COMPOUND_CLASS_ANNOTATION = "compoundClassAnnotation";
  private CompoundClasses compoundClassAnnotation;

  public static final String JSON_PROPERTY_CONFIDENCE_EXACT_MATCH = "confidenceExactMatch";
  private Double confidenceExactMatch;

  public static final String JSON_PROPERTY_CONFIDENCE_APPROX_MATCH = "confidenceApproxMatch";
  private Double confidenceApproxMatch;

  public static final String JSON_PROPERTY_EXPANSIVE_SEARCH_STATE = "expansiveSearchState";
  private ConfidenceMode expansiveSearchState;

  public static final String JSON_PROPERTY_SPECIFIED_DATABASES = "specifiedDatabases";
  private List<String> specifiedDatabases;

  public static final String JSON_PROPERTY_EXPANDED_DATABASES = "expandedDatabases";
  private List<String> expandedDatabases;

  public FeatureAnnotations() {
  }

  public FeatureAnnotations formulaAnnotation(FormulaCandidate formulaAnnotation) {
    
    this.formulaAnnotation = formulaAnnotation;
    return this;
  }

   /**
   * Get formulaAnnotation
   * @return formulaAnnotation
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_FORMULA_ANNOTATION)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public FormulaCandidate getFormulaAnnotation() {
    return formulaAnnotation;
  }


  @JsonProperty(JSON_PROPERTY_FORMULA_ANNOTATION)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setFormulaAnnotation(FormulaCandidate formulaAnnotation) {
    this.formulaAnnotation = formulaAnnotation;
  }

  public FeatureAnnotations structureAnnotation(StructureCandidateScored structureAnnotation) {
    
    this.structureAnnotation = structureAnnotation;
    return this;
  }

   /**
   * Get structureAnnotation
   * @return structureAnnotation
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_STRUCTURE_ANNOTATION)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public StructureCandidateScored getStructureAnnotation() {
    return structureAnnotation;
  }


  @JsonProperty(JSON_PROPERTY_STRUCTURE_ANNOTATION)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setStructureAnnotation(StructureCandidateScored structureAnnotation) {
    this.structureAnnotation = structureAnnotation;
  }

  public FeatureAnnotations compoundClassAnnotation(CompoundClasses compoundClassAnnotation) {
    
    this.compoundClassAnnotation = compoundClassAnnotation;
    return this;
  }

   /**
   * Get compoundClassAnnotation
   * @return compoundClassAnnotation
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_COMPOUND_CLASS_ANNOTATION)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public CompoundClasses getCompoundClassAnnotation() {
    return compoundClassAnnotation;
  }


  @JsonProperty(JSON_PROPERTY_COMPOUND_CLASS_ANNOTATION)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setCompoundClassAnnotation(CompoundClasses compoundClassAnnotation) {
    this.compoundClassAnnotation = compoundClassAnnotation;
  }

  public FeatureAnnotations confidenceExactMatch(Double confidenceExactMatch) {
    
    this.confidenceExactMatch = confidenceExactMatch;
    return this;
  }

   /**
   * Confidence Score that represents the confidence whether the top hit is correct.
   * @return confidenceExactMatch
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_CONFIDENCE_EXACT_MATCH)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public Double getConfidenceExactMatch() {
    return confidenceExactMatch;
  }


  @JsonProperty(JSON_PROPERTY_CONFIDENCE_EXACT_MATCH)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setConfidenceExactMatch(Double confidenceExactMatch) {
    this.confidenceExactMatch = confidenceExactMatch;
  }

  public FeatureAnnotations confidenceApproxMatch(Double confidenceApproxMatch) {
    
    this.confidenceApproxMatch = confidenceApproxMatch;
    return this;
  }

   /**
   * Confidence Score that represents the confidence whether the top hit or a very similar hit (estimated by MCES distance) is correct.
   * @return confidenceApproxMatch
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_CONFIDENCE_APPROX_MATCH)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public Double getConfidenceApproxMatch() {
    return confidenceApproxMatch;
  }


  @JsonProperty(JSON_PROPERTY_CONFIDENCE_APPROX_MATCH)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setConfidenceApproxMatch(Double confidenceApproxMatch) {
    this.confidenceApproxMatch = confidenceApproxMatch;
  }

  public FeatureAnnotations expansiveSearchState(ConfidenceMode expansiveSearchState) {
    
    this.expansiveSearchState = expansiveSearchState;
    return this;
  }

   /**
   * Get expansiveSearchState
   * @return expansiveSearchState
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_EXPANSIVE_SEARCH_STATE)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public ConfidenceMode getExpansiveSearchState() {
    return expansiveSearchState;
  }


  @JsonProperty(JSON_PROPERTY_EXPANSIVE_SEARCH_STATE)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setExpansiveSearchState(ConfidenceMode expansiveSearchState) {
    this.expansiveSearchState = expansiveSearchState;
  }

  public FeatureAnnotations specifiedDatabases(List<String> specifiedDatabases) {
    
    this.specifiedDatabases = specifiedDatabases;
    return this;
  }

  public FeatureAnnotations addSpecifiedDatabasesItem(String specifiedDatabasesItem) {
    if (this.specifiedDatabases == null) {
      this.specifiedDatabases = new ArrayList<>();
    }
    this.specifiedDatabases.add(specifiedDatabasesItem);
    return this;
  }

   /**
   * List of databases that have been specified by for structure db search. Null if no structure db search has been performed.
   * @return specifiedDatabases
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_SPECIFIED_DATABASES)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public List<String> getSpecifiedDatabases() {
    return specifiedDatabases;
  }


  @JsonProperty(JSON_PROPERTY_SPECIFIED_DATABASES)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setSpecifiedDatabases(List<String> specifiedDatabases) {
    this.specifiedDatabases = specifiedDatabases;
  }

  public FeatureAnnotations expandedDatabases(List<String> expandedDatabases) {
    
    this.expandedDatabases = expandedDatabases;
    return this;
  }

  public FeatureAnnotations addExpandedDatabasesItem(String expandedDatabasesItem) {
    if (this.expandedDatabases == null) {
      this.expandedDatabases = new ArrayList<>();
    }
    this.expandedDatabases.add(expandedDatabasesItem);
    return this;
  }

   /**
   * List of databases that have been used to expand search space during expansive search. Null if no structure db search has been performed.
   * @return expandedDatabases
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_EXPANDED_DATABASES)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public List<String> getExpandedDatabases() {
    return expandedDatabases;
  }


  @JsonProperty(JSON_PROPERTY_EXPANDED_DATABASES)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setExpandedDatabases(List<String> expandedDatabases) {
    this.expandedDatabases = expandedDatabases;
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    FeatureAnnotations featureAnnotations = (FeatureAnnotations) o;
    return Objects.equals(this.formulaAnnotation, featureAnnotations.formulaAnnotation) &&
        Objects.equals(this.structureAnnotation, featureAnnotations.structureAnnotation) &&
        Objects.equals(this.compoundClassAnnotation, featureAnnotations.compoundClassAnnotation) &&
        Objects.equals(this.confidenceExactMatch, featureAnnotations.confidenceExactMatch) &&
        Objects.equals(this.confidenceApproxMatch, featureAnnotations.confidenceApproxMatch) &&
        Objects.equals(this.expansiveSearchState, featureAnnotations.expansiveSearchState) &&
        Objects.equals(this.specifiedDatabases, featureAnnotations.specifiedDatabases) &&
        Objects.equals(this.expandedDatabases, featureAnnotations.expandedDatabases);
  }

  @Override
  public int hashCode() {
    return Objects.hash(formulaAnnotation, structureAnnotation, compoundClassAnnotation, confidenceExactMatch, confidenceApproxMatch, expansiveSearchState, specifiedDatabases, expandedDatabases);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class FeatureAnnotations {\n");
    sb.append("    formulaAnnotation: ").append(toIndentedString(formulaAnnotation)).append("\n");
    sb.append("    structureAnnotation: ").append(toIndentedString(structureAnnotation)).append("\n");
    sb.append("    compoundClassAnnotation: ").append(toIndentedString(compoundClassAnnotation)).append("\n");
    sb.append("    confidenceExactMatch: ").append(toIndentedString(confidenceExactMatch)).append("\n");
    sb.append("    confidenceApproxMatch: ").append(toIndentedString(confidenceApproxMatch)).append("\n");
    sb.append("    expansiveSearchState: ").append(toIndentedString(expansiveSearchState)).append("\n");
    sb.append("    specifiedDatabases: ").append(toIndentedString(specifiedDatabases)).append("\n");
    sb.append("    expandedDatabases: ").append(toIndentedString(expandedDatabases)).append("\n");
    sb.append("}");
    return sb.toString();
  }

  /**
   * Convert the given object to string with each line indented by 4 spaces
   * (except the first line).
   */
  private String toIndentedString(Object o) {
    if (o == null) {
      return "null";
    }
    return o.toString().replace("\n", "\n    ");
  }

}

