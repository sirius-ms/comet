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
import com.fasterxml.jackson.annotation.JsonPropertyOrder;
import com.fasterxml.jackson.annotation.JsonTypeName;

/**
 * 
 */
@JsonPropertyOrder({
  LossEdge.JSON_PROPERTY_SOURCE_FRAGMENT,
  LossEdge.JSON_PROPERTY_TARGET_FRAGMENT,
  LossEdge.JSON_PROPERTY_MOLECULAR_FORMULA,
  LossEdge.JSON_PROPERTY_SCORE
})
@javax.annotation.Generated(value = "org.openapitools.codegen.languages.JavaClientCodegen")
public class LossEdge {
  public static final String JSON_PROPERTY_SOURCE_FRAGMENT = "sourceFragment";
  private FragmentNode sourceFragment;

  public static final String JSON_PROPERTY_TARGET_FRAGMENT = "targetFragment";
  private FragmentNode targetFragment;

  public static final String JSON_PROPERTY_MOLECULAR_FORMULA = "molecularFormula";
  private String molecularFormula;

  public static final String JSON_PROPERTY_SCORE = "score";
  private Double score;

  public LossEdge() {
  }

  public LossEdge sourceFragment(FragmentNode sourceFragment) {
    
    this.sourceFragment = sourceFragment;
    return this;
  }

   /**
   * Get sourceFragment
   * @return sourceFragment
  **/
  @javax.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_SOURCE_FRAGMENT)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public FragmentNode getSourceFragment() {
    return sourceFragment;
  }


  @JsonProperty(JSON_PROPERTY_SOURCE_FRAGMENT)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setSourceFragment(FragmentNode sourceFragment) {
    this.sourceFragment = sourceFragment;
  }


  public LossEdge targetFragment(FragmentNode targetFragment) {
    
    this.targetFragment = targetFragment;
    return this;
  }

   /**
   * Get targetFragment
   * @return targetFragment
  **/
  @javax.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_TARGET_FRAGMENT)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public FragmentNode getTargetFragment() {
    return targetFragment;
  }


  @JsonProperty(JSON_PROPERTY_TARGET_FRAGMENT)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setTargetFragment(FragmentNode targetFragment) {
    this.targetFragment = targetFragment;
  }


  public LossEdge molecularFormula(String molecularFormula) {
    
    this.molecularFormula = molecularFormula;
    return this;
  }

   /**
   * Get molecularFormula
   * @return molecularFormula
  **/
  @javax.annotation.Nullable
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


  public LossEdge score(Double score) {
    
    this.score = score;
    return this;
  }

   /**
   * Get score
   * @return score
  **/
  @javax.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_SCORE)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public Double getScore() {
    return score;
  }


  @JsonProperty(JSON_PROPERTY_SCORE)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setScore(Double score) {
    this.score = score;
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    LossEdge lossEdge = (LossEdge) o;
    return Objects.equals(this.sourceFragment, lossEdge.sourceFragment) &&
        Objects.equals(this.targetFragment, lossEdge.targetFragment) &&
        Objects.equals(this.molecularFormula, lossEdge.molecularFormula) &&
        Objects.equals(this.score, lossEdge.score);
  }

  @Override
  public int hashCode() {
    return Objects.hash(sourceFragment, targetFragment, molecularFormula, score);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class LossEdge {\n");
    sb.append("    sourceFragment: ").append(toIndentedString(sourceFragment)).append("\n");
    sb.append("    targetFragment: ").append(toIndentedString(targetFragment)).append("\n");
    sb.append("    molecularFormula: ").append(toIndentedString(molecularFormula)).append("\n");
    sb.append("    score: ").append(toIndentedString(score)).append("\n");
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

