/*
 * SIRIUS Nightsky API
 * REST API that provides the full functionality of SIRIUS and its web services as background service. It is intended as entry-point for scripting languages and software integration SDKs.This API is exposed by SIRIUS 6.0.0-SNAPSHOT
 *
 * The version of the OpenAPI document: 2.0
 * 
 *
 * NOTE: This class is auto generated by OpenAPI Generator (https://openapi-generator.tech).
 * https://openapi-generator.tech
 * Do not edit the class manually.
 */


package de.unijena.bioinf.ms.nightsky.sdk.model;

import java.util.Objects;
import java.util.Arrays;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonTypeName;
import com.fasterxml.jackson.annotation.JsonValue;
import de.unijena.bioinf.ms.nightsky.sdk.model.FragmentNode;
import de.unijena.bioinf.ms.nightsky.sdk.model.LossEdge;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;
import com.fasterxml.jackson.annotation.JsonTypeName;

/**
 * 
 */
@JsonPropertyOrder({
  FragmentationTree.JSON_PROPERTY_FRAGMENTS,
  FragmentationTree.JSON_PROPERTY_LOSSES,
  FragmentationTree.JSON_PROPERTY_TREE_SCORE,
  FragmentationTree.JSON_PROPERTY_ROOT
})
@javax.annotation.Generated(value = "org.openapitools.codegen.languages.JavaClientCodegen")
public class FragmentationTree {
  public static final String JSON_PROPERTY_FRAGMENTS = "fragments";
  private List<FragmentNode> fragments;

  public static final String JSON_PROPERTY_LOSSES = "losses";
  private List<LossEdge> losses;

  public static final String JSON_PROPERTY_TREE_SCORE = "treeScore";
  private Double treeScore;

  public static final String JSON_PROPERTY_ROOT = "root";
  private FragmentNode root;

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
  @javax.annotation.Nullable
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
  @javax.annotation.Nullable
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
  @javax.annotation.Nullable
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


  public FragmentationTree root(FragmentNode root) {
    
    this.root = root;
    return this;
  }

   /**
   * Get root
   * @return root
  **/
  @javax.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_ROOT)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public FragmentNode getRoot() {
    return root;
  }


  @JsonProperty(JSON_PROPERTY_ROOT)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setRoot(FragmentNode root) {
    this.root = root;
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
        Objects.equals(this.root, fragmentationTree.root);
  }

  @Override
  public int hashCode() {
    return Objects.hash(fragments, losses, treeScore, root);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class FragmentationTree {\n");
    sb.append("    fragments: ").append(toIndentedString(fragments)).append("\n");
    sb.append("    losses: ").append(toIndentedString(losses)).append("\n");
    sb.append("    treeScore: ").append(toIndentedString(treeScore)).append("\n");
    sb.append("    root: ").append(toIndentedString(root)).append("\n");
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

