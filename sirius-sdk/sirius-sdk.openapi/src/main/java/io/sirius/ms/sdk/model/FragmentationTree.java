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
import io.sirius.ms.sdk.model.FragmentNode;
import io.sirius.ms.sdk.model.LossEdge;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;
import com.fasterxml.jackson.annotation.JsonTypeName;

/**
 * Simple and easy serializable fragmentation tree model with annotated fragments/nodes abd losses/edges  Root fragment has index 0;  Molecular formula and adduct are identical to the ones of the corresponding molecular formula candidate and SpectrumAnnotation
 */
@JsonPropertyOrder({
  FragmentationTree.JSON_PROPERTY_FRAGMENTS,
  FragmentationTree.JSON_PROPERTY_LOSSES,
  FragmentationTree.JSON_PROPERTY_TREE_SCORE,
  FragmentationTree.JSON_PROPERTY_MOLECULAR_FORMULA,
  FragmentationTree.JSON_PROPERTY_ADDUCT
})
@jakarta.annotation.Generated(value = "org.openapitools.codegen.languages.JavaClientCodegen", comments = "Generator version: 7.6.0")
public class FragmentationTree {
  public static final String JSON_PROPERTY_FRAGMENTS = "fragments";
  private List<FragmentNode> fragments = new ArrayList<>();

  public static final String JSON_PROPERTY_LOSSES = "losses";
  private List<LossEdge> losses = new ArrayList<>();

  public static final String JSON_PROPERTY_TREE_SCORE = "treeScore";
  private Double treeScore;

  public static final String JSON_PROPERTY_MOLECULAR_FORMULA = "molecularFormula";
  private String molecularFormula;

  public static final String JSON_PROPERTY_ADDUCT = "adduct";
  private String adduct;

  public FragmentationTree() {
  }

  public FragmentationTree fragments(List<FragmentNode> fragments) {
    
    this.fragments = fragments;
    return this;
  }

  public FragmentationTree addFragmentsItem(FragmentNode fragmentsItem) {
    if (this.fragments == null) {
      this.fragments = new ArrayList<>();
    }
    this.fragments.add(fragmentsItem);
    return this;
  }

   /**
   * Get fragments
   * @return fragments
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_FRAGMENTS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public List<FragmentNode> getFragments() {
    return fragments;
  }


  @JsonProperty(JSON_PROPERTY_FRAGMENTS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setFragments(List<FragmentNode> fragments) {
    this.fragments = fragments;
  }

  public FragmentationTree losses(List<LossEdge> losses) {
    
    this.losses = losses;
    return this;
  }

  public FragmentationTree addLossesItem(LossEdge lossesItem) {
    if (this.losses == null) {
      this.losses = new ArrayList<>();
    }
    this.losses.add(lossesItem);
    return this;
  }

   /**
   * Get losses
   * @return losses
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_LOSSES)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public List<LossEdge> getLosses() {
    return losses;
  }


  @JsonProperty(JSON_PROPERTY_LOSSES)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setLosses(List<LossEdge> losses) {
    this.losses = losses;
  }

  public FragmentationTree treeScore(Double treeScore) {
    
    this.treeScore = treeScore;
    return this;
  }

   /**
   * Get treeScore
   * @return treeScore
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_TREE_SCORE)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public Double getTreeScore() {
    return treeScore;
  }


  @JsonProperty(JSON_PROPERTY_TREE_SCORE)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setTreeScore(Double treeScore) {
    this.treeScore = treeScore;
  }

  public FragmentationTree molecularFormula(String molecularFormula) {
    
    this.molecularFormula = molecularFormula;
    return this;
  }

   /**
   * Get molecularFormula
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

  public FragmentationTree adduct(String adduct) {
    
    this.adduct = adduct;
    return this;
  }

   /**
   * Get adduct
   * @return adduct
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_ADDUCT)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public String getAdduct() {
    return adduct;
  }


  @JsonProperty(JSON_PROPERTY_ADDUCT)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setAdduct(String adduct) {
    this.adduct = adduct;
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    FragmentationTree fragmentationTree = (FragmentationTree) o;
    return Objects.equals(this.fragments, fragmentationTree.fragments) &&
        Objects.equals(this.losses, fragmentationTree.losses) &&
        Objects.equals(this.treeScore, fragmentationTree.treeScore) &&
        Objects.equals(this.molecularFormula, fragmentationTree.molecularFormula) &&
        Objects.equals(this.adduct, fragmentationTree.adduct);
  }

  @Override
  public int hashCode() {
    return Objects.hash(fragments, losses, treeScore, molecularFormula, adduct);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class FragmentationTree {\n");
    sb.append("    fragments: ").append(toIndentedString(fragments)).append("\n");
    sb.append("    losses: ").append(toIndentedString(losses)).append("\n");
    sb.append("    treeScore: ").append(toIndentedString(treeScore)).append("\n");
    sb.append("    molecularFormula: ").append(toIndentedString(molecularFormula)).append("\n");
    sb.append("    adduct: ").append(toIndentedString(adduct)).append("\n");
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

