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
import io.sirius.ms.sdk.model.ConsensusCriterionCSI;
import io.sirius.ms.sdk.model.StructureCandidate;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;
import com.fasterxml.jackson.annotation.JsonTypeName;

/**
 * ConsensusAnnotationsCSI
 */
@JsonPropertyOrder({
  ConsensusAnnotationsCSI.JSON_PROPERTY_MOLECULAR_FORMULA,
  ConsensusAnnotationsCSI.JSON_PROPERTY_COMPOUND_CLASSES,
  ConsensusAnnotationsCSI.JSON_PROPERTY_SUPPORTING_FEATURE_IDS,
  ConsensusAnnotationsCSI.JSON_PROPERTY_SELECTION_CRITERION,
  ConsensusAnnotationsCSI.JSON_PROPERTY_CSI_FINGER_ID_STRUCTURE,
  ConsensusAnnotationsCSI.JSON_PROPERTY_CONFIDENCE_EXACT_MATCH,
  ConsensusAnnotationsCSI.JSON_PROPERTY_CONFIDENCE_APPROX_MATCH
})
@jakarta.annotation.Generated(value = "org.openapitools.codegen.languages.JavaClientCodegen", comments = "Generator version: 7.6.0")
public class ConsensusAnnotationsCSI {
  public static final String JSON_PROPERTY_MOLECULAR_FORMULA = "molecularFormula";
  private String molecularFormula;

  public static final String JSON_PROPERTY_COMPOUND_CLASSES = "compoundClasses";
  private CompoundClasses compoundClasses;

  public static final String JSON_PROPERTY_SUPPORTING_FEATURE_IDS = "supportingFeatureIds";
  private List<String> supportingFeatureIds;

  public static final String JSON_PROPERTY_SELECTION_CRITERION = "selectionCriterion";
  private ConsensusCriterionCSI selectionCriterion;

  public static final String JSON_PROPERTY_CSI_FINGER_ID_STRUCTURE = "csiFingerIdStructure";
  private StructureCandidate csiFingerIdStructure;

  public static final String JSON_PROPERTY_CONFIDENCE_EXACT_MATCH = "confidenceExactMatch";
  private Double confidenceExactMatch;

  public static final String JSON_PROPERTY_CONFIDENCE_APPROX_MATCH = "confidenceApproxMatch";
  private Double confidenceApproxMatch;

  public ConsensusAnnotationsCSI() {
  }

  public ConsensusAnnotationsCSI molecularFormula(String molecularFormula) {
    
    this.molecularFormula = molecularFormula;
    return this;
  }

   /**
   * Molecular formula of the consensus annotation  Might be null if no consensus formula is available.
   * @return molecularFormula
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_MOLECULAR_FORMULA)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public String getMolecularFormula() {
    return molecularFormula;
  }


  @JsonProperty(JSON_PROPERTY_MOLECULAR_FORMULA)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setMolecularFormula(String molecularFormula) {
    this.molecularFormula = molecularFormula;
  }

  public ConsensusAnnotationsCSI compoundClasses(CompoundClasses compoundClasses) {
    
    this.compoundClasses = compoundClasses;
    return this;
  }

   /**
   * Get compoundClasses
   * @return compoundClasses
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_COMPOUND_CLASSES)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public CompoundClasses getCompoundClasses() {
    return compoundClasses;
  }


  @JsonProperty(JSON_PROPERTY_COMPOUND_CLASSES)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setCompoundClasses(CompoundClasses compoundClasses) {
    this.compoundClasses = compoundClasses;
  }

  public ConsensusAnnotationsCSI supportingFeatureIds(List<String> supportingFeatureIds) {
    
    this.supportingFeatureIds = supportingFeatureIds;
    return this;
  }

  public ConsensusAnnotationsCSI addSupportingFeatureIdsItem(String supportingFeatureIdsItem) {
    if (this.supportingFeatureIds == null) {
      this.supportingFeatureIds = new ArrayList<>();
    }
    this.supportingFeatureIds.add(supportingFeatureIdsItem);
    return this;
  }

   /**
   * FeatureIds where the topAnnotation supports this annotation.
   * @return supportingFeatureIds
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_SUPPORTING_FEATURE_IDS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public List<String> getSupportingFeatureIds() {
    return supportingFeatureIds;
  }


  @JsonProperty(JSON_PROPERTY_SUPPORTING_FEATURE_IDS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setSupportingFeatureIds(List<String> supportingFeatureIds) {
    this.supportingFeatureIds = supportingFeatureIds;
  }

  public ConsensusAnnotationsCSI selectionCriterion(ConsensusCriterionCSI selectionCriterion) {
    
    this.selectionCriterion = selectionCriterion;
    return this;
  }

   /**
   * Get selectionCriterion
   * @return selectionCriterion
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_SELECTION_CRITERION)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public ConsensusCriterionCSI getSelectionCriterion() {
    return selectionCriterion;
  }


  @JsonProperty(JSON_PROPERTY_SELECTION_CRITERION)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setSelectionCriterion(ConsensusCriterionCSI selectionCriterion) {
    this.selectionCriterion = selectionCriterion;
  }

  public ConsensusAnnotationsCSI csiFingerIdStructure(StructureCandidate csiFingerIdStructure) {
    
    this.csiFingerIdStructure = csiFingerIdStructure;
    return this;
  }

   /**
   * Get csiFingerIdStructure
   * @return csiFingerIdStructure
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_CSI_FINGER_ID_STRUCTURE)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public StructureCandidate getCsiFingerIdStructure() {
    return csiFingerIdStructure;
  }


  @JsonProperty(JSON_PROPERTY_CSI_FINGER_ID_STRUCTURE)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setCsiFingerIdStructure(StructureCandidate csiFingerIdStructure) {
    this.csiFingerIdStructure = csiFingerIdStructure;
  }

  public ConsensusAnnotationsCSI confidenceExactMatch(Double confidenceExactMatch) {
    
    this.confidenceExactMatch = confidenceExactMatch;
    return this;
  }

   /**
   * Confidence value that represents the certainty that reported consensus structure is exactly the measured one  If multiple features support this consensus structure the maximum confidence is reported
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

  public ConsensusAnnotationsCSI confidenceApproxMatch(Double confidenceApproxMatch) {
    
    this.confidenceApproxMatch = confidenceApproxMatch;
    return this;
  }

   /**
   * Confidence value that represents the certainty that the exact consensus structure or a very similar  structure (e.g. measured by Maximum Common Edge Subgraph Distance) is the measured one.  If multiple features support this consensus structure the maximum confidence is reported
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

  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    ConsensusAnnotationsCSI consensusAnnotationsCSI = (ConsensusAnnotationsCSI) o;
    return Objects.equals(this.molecularFormula, consensusAnnotationsCSI.molecularFormula) &&
        Objects.equals(this.compoundClasses, consensusAnnotationsCSI.compoundClasses) &&
        Objects.equals(this.supportingFeatureIds, consensusAnnotationsCSI.supportingFeatureIds) &&
        Objects.equals(this.selectionCriterion, consensusAnnotationsCSI.selectionCriterion) &&
        Objects.equals(this.csiFingerIdStructure, consensusAnnotationsCSI.csiFingerIdStructure) &&
        Objects.equals(this.confidenceExactMatch, consensusAnnotationsCSI.confidenceExactMatch) &&
        Objects.equals(this.confidenceApproxMatch, consensusAnnotationsCSI.confidenceApproxMatch);
  }

  @Override
  public int hashCode() {
    return Objects.hash(molecularFormula, compoundClasses, supportingFeatureIds, selectionCriterion, csiFingerIdStructure, confidenceExactMatch, confidenceApproxMatch);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class ConsensusAnnotationsCSI {\n");
    sb.append("    molecularFormula: ").append(toIndentedString(molecularFormula)).append("\n");
    sb.append("    compoundClasses: ").append(toIndentedString(compoundClasses)).append("\n");
    sb.append("    supportingFeatureIds: ").append(toIndentedString(supportingFeatureIds)).append("\n");
    sb.append("    selectionCriterion: ").append(toIndentedString(selectionCriterion)).append("\n");
    sb.append("    csiFingerIdStructure: ").append(toIndentedString(csiFingerIdStructure)).append("\n");
    sb.append("    confidenceExactMatch: ").append(toIndentedString(confidenceExactMatch)).append("\n");
    sb.append("    confidenceApproxMatch: ").append(toIndentedString(confidenceApproxMatch)).append("\n");
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

