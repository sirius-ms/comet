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
import io.sirius.ms.sdk.model.ZodiacEdgeFilterThresholds;
import io.sirius.ms.sdk.model.ZodiacEpochs;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;
import com.fasterxml.jackson.annotation.JsonTypeName;

/**
 * User/developer friendly parameter subset for the ZODIAC tool (Network base molecular formula re-ranking).  Needs results from Formula/SIRIUS Tool
 */
@JsonPropertyOrder({
  Zodiac.JSON_PROPERTY_ENABLED,
  Zodiac.JSON_PROPERTY_CONSIDERED_CANDIDATES_AT300_MZ,
  Zodiac.JSON_PROPERTY_CONSIDERED_CANDIDATES_AT800_MZ,
  Zodiac.JSON_PROPERTY_RUN_IN_TWO_STEPS,
  Zodiac.JSON_PROPERTY_EDGE_FILTER_THRESHOLDS,
  Zodiac.JSON_PROPERTY_GIBBS_SAMPLER_PARAMETERS
})
@jakarta.annotation.Generated(value = "org.openapitools.codegen.languages.JavaClientCodegen", comments = "Generator version: 7.6.0")
public class Zodiac {
  public static final String JSON_PROPERTY_ENABLED = "enabled";
  private Boolean enabled;

  public static final String JSON_PROPERTY_CONSIDERED_CANDIDATES_AT300_MZ = "consideredCandidatesAt300Mz";
  private Integer consideredCandidatesAt300Mz;

  public static final String JSON_PROPERTY_CONSIDERED_CANDIDATES_AT800_MZ = "consideredCandidatesAt800Mz";
  private Integer consideredCandidatesAt800Mz;

  public static final String JSON_PROPERTY_RUN_IN_TWO_STEPS = "runInTwoSteps";
  private Boolean runInTwoSteps;

  public static final String JSON_PROPERTY_EDGE_FILTER_THRESHOLDS = "edgeFilterThresholds";
  private ZodiacEdgeFilterThresholds edgeFilterThresholds;

  public static final String JSON_PROPERTY_GIBBS_SAMPLER_PARAMETERS = "gibbsSamplerParameters";
  private ZodiacEpochs gibbsSamplerParameters;

  public Zodiac() {
  }

  public Zodiac enabled(Boolean enabled) {
    
    this.enabled = enabled;
    return this;
  }

   /**
   * tags whether the tool is enabled
   * @return enabled
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_ENABLED)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public Boolean isEnabled() {
    return enabled;
  }


  @JsonProperty(JSON_PROPERTY_ENABLED)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setEnabled(Boolean enabled) {
    this.enabled = enabled;
  }

  public Zodiac consideredCandidatesAt300Mz(Integer consideredCandidatesAt300Mz) {
    
    this.consideredCandidatesAt300Mz = consideredCandidatesAt300Mz;
    return this;
  }

   /**
   * Maximum number of candidate molecular formulas (fragmentation trees computed by SIRIUS) per compound which are considered by ZODIAC for compounds below 300 m/z.
   * @return consideredCandidatesAt300Mz
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_CONSIDERED_CANDIDATES_AT300_MZ)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public Integer getConsideredCandidatesAt300Mz() {
    return consideredCandidatesAt300Mz;
  }


  @JsonProperty(JSON_PROPERTY_CONSIDERED_CANDIDATES_AT300_MZ)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setConsideredCandidatesAt300Mz(Integer consideredCandidatesAt300Mz) {
    this.consideredCandidatesAt300Mz = consideredCandidatesAt300Mz;
  }

  public Zodiac consideredCandidatesAt800Mz(Integer consideredCandidatesAt800Mz) {
    
    this.consideredCandidatesAt800Mz = consideredCandidatesAt800Mz;
    return this;
  }

   /**
   * Maximum number of candidate molecular formulas (fragmentation trees computed by SIRIUS) per compound which are considered by ZODIAC for compounds above 800 m/z.
   * @return consideredCandidatesAt800Mz
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_CONSIDERED_CANDIDATES_AT800_MZ)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public Integer getConsideredCandidatesAt800Mz() {
    return consideredCandidatesAt800Mz;
  }


  @JsonProperty(JSON_PROPERTY_CONSIDERED_CANDIDATES_AT800_MZ)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setConsideredCandidatesAt800Mz(Integer consideredCandidatesAt800Mz) {
    this.consideredCandidatesAt800Mz = consideredCandidatesAt800Mz;
  }

  public Zodiac runInTwoSteps(Boolean runInTwoSteps) {
    
    this.runInTwoSteps = runInTwoSteps;
    return this;
  }

   /**
   * As default ZODIAC runs a 2-step approach. First running &#39;good quality compounds&#39; only, and afterwards including the remaining.
   * @return runInTwoSteps
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_RUN_IN_TWO_STEPS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public Boolean isRunInTwoSteps() {
    return runInTwoSteps;
  }


  @JsonProperty(JSON_PROPERTY_RUN_IN_TWO_STEPS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setRunInTwoSteps(Boolean runInTwoSteps) {
    this.runInTwoSteps = runInTwoSteps;
  }

  public Zodiac edgeFilterThresholds(ZodiacEdgeFilterThresholds edgeFilterThresholds) {
    
    this.edgeFilterThresholds = edgeFilterThresholds;
    return this;
  }

   /**
   * Get edgeFilterThresholds
   * @return edgeFilterThresholds
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_EDGE_FILTER_THRESHOLDS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public ZodiacEdgeFilterThresholds getEdgeFilterThresholds() {
    return edgeFilterThresholds;
  }


  @JsonProperty(JSON_PROPERTY_EDGE_FILTER_THRESHOLDS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setEdgeFilterThresholds(ZodiacEdgeFilterThresholds edgeFilterThresholds) {
    this.edgeFilterThresholds = edgeFilterThresholds;
  }

  public Zodiac gibbsSamplerParameters(ZodiacEpochs gibbsSamplerParameters) {
    
    this.gibbsSamplerParameters = gibbsSamplerParameters;
    return this;
  }

   /**
   * Get gibbsSamplerParameters
   * @return gibbsSamplerParameters
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_GIBBS_SAMPLER_PARAMETERS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public ZodiacEpochs getGibbsSamplerParameters() {
    return gibbsSamplerParameters;
  }


  @JsonProperty(JSON_PROPERTY_GIBBS_SAMPLER_PARAMETERS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setGibbsSamplerParameters(ZodiacEpochs gibbsSamplerParameters) {
    this.gibbsSamplerParameters = gibbsSamplerParameters;
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    Zodiac zodiac = (Zodiac) o;
    return Objects.equals(this.enabled, zodiac.enabled) &&
        Objects.equals(this.consideredCandidatesAt300Mz, zodiac.consideredCandidatesAt300Mz) &&
        Objects.equals(this.consideredCandidatesAt800Mz, zodiac.consideredCandidatesAt800Mz) &&
        Objects.equals(this.runInTwoSteps, zodiac.runInTwoSteps) &&
        Objects.equals(this.edgeFilterThresholds, zodiac.edgeFilterThresholds) &&
        Objects.equals(this.gibbsSamplerParameters, zodiac.gibbsSamplerParameters);
  }

  @Override
  public int hashCode() {
    return Objects.hash(enabled, consideredCandidatesAt300Mz, consideredCandidatesAt800Mz, runInTwoSteps, edgeFilterThresholds, gibbsSamplerParameters);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class Zodiac {\n");
    sb.append("    enabled: ").append(toIndentedString(enabled)).append("\n");
    sb.append("    consideredCandidatesAt300Mz: ").append(toIndentedString(consideredCandidatesAt300Mz)).append("\n");
    sb.append("    consideredCandidatesAt800Mz: ").append(toIndentedString(consideredCandidatesAt800Mz)).append("\n");
    sb.append("    runInTwoSteps: ").append(toIndentedString(runInTwoSteps)).append("\n");
    sb.append("    edgeFilterThresholds: ").append(toIndentedString(edgeFilterThresholds)).append("\n");
    sb.append("    gibbsSamplerParameters: ").append(toIndentedString(gibbsSamplerParameters)).append("\n");
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

