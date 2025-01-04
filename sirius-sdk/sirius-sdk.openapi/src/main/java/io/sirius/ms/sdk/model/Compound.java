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
import io.sirius.ms.sdk.model.AlignedFeature;
import io.sirius.ms.sdk.model.ConsensusAnnotationsCSI;
import io.sirius.ms.sdk.model.ConsensusAnnotationsDeNovo;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;
import com.fasterxml.jackson.annotation.JsonTypeName;

/**
 * Compound
 */
@JsonPropertyOrder({
  Compound.JSON_PROPERTY_COMPOUND_ID,
  Compound.JSON_PROPERTY_NAME,
  Compound.JSON_PROPERTY_RT_START_SECONDS,
  Compound.JSON_PROPERTY_RT_END_SECONDS,
  Compound.JSON_PROPERTY_NEUTRAL_MASS,
  Compound.JSON_PROPERTY_FEATURES,
  Compound.JSON_PROPERTY_CONSENSUS_ANNOTATIONS,
  Compound.JSON_PROPERTY_CONSENSUS_ANNOTATIONS_DE_NOVO,
  Compound.JSON_PROPERTY_CUSTOM_ANNOTATIONS
})
@jakarta.annotation.Generated(value = "org.openapitools.codegen.languages.JavaClientCodegen", comments = "Generator version: 7.6.0")
public class Compound {
  public static final String JSON_PROPERTY_COMPOUND_ID = "compoundId";
  private String compoundId;

  public static final String JSON_PROPERTY_NAME = "name";
  private String name;

  public static final String JSON_PROPERTY_RT_START_SECONDS = "rtStartSeconds";
  private Double rtStartSeconds;

  public static final String JSON_PROPERTY_RT_END_SECONDS = "rtEndSeconds";
  private Double rtEndSeconds;

  public static final String JSON_PROPERTY_NEUTRAL_MASS = "neutralMass";
  private Double neutralMass;

  public static final String JSON_PROPERTY_FEATURES = "features";
  private List<AlignedFeature> features = new ArrayList<>();

  public static final String JSON_PROPERTY_CONSENSUS_ANNOTATIONS = "consensusAnnotations";
  private ConsensusAnnotationsCSI consensusAnnotations;

  public static final String JSON_PROPERTY_CONSENSUS_ANNOTATIONS_DE_NOVO = "consensusAnnotationsDeNovo";
  private ConsensusAnnotationsDeNovo consensusAnnotationsDeNovo;

  public static final String JSON_PROPERTY_CUSTOM_ANNOTATIONS = "customAnnotations";
  private ConsensusAnnotationsCSI customAnnotations;

  public Compound() {
  }

  public Compound compoundId(String compoundId) {
    
    this.compoundId = compoundId;
    return this;
  }

   /**
   * uid of this compound Entity
   * @return compoundId
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_COMPOUND_ID)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public String getCompoundId() {
    return compoundId;
  }


  @JsonProperty(JSON_PROPERTY_COMPOUND_ID)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setCompoundId(String compoundId) {
    this.compoundId = compoundId;
  }

  public Compound name(String name) {
    
    this.name = name;
    return this;
  }

   /**
   * Some (optional) human-readable name
   * @return name
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_NAME)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public String getName() {
    return name;
  }


  @JsonProperty(JSON_PROPERTY_NAME)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setName(String name) {
    this.name = name;
  }

  public Compound rtStartSeconds(Double rtStartSeconds) {
    
    this.rtStartSeconds = rtStartSeconds;
    return this;
  }

   /**
   * The merged/consensus retention time start (earliest rt) of this compound
   * @return rtStartSeconds
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_RT_START_SECONDS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public Double getRtStartSeconds() {
    return rtStartSeconds;
  }


  @JsonProperty(JSON_PROPERTY_RT_START_SECONDS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setRtStartSeconds(Double rtStartSeconds) {
    this.rtStartSeconds = rtStartSeconds;
  }

  public Compound rtEndSeconds(Double rtEndSeconds) {
    
    this.rtEndSeconds = rtEndSeconds;
    return this;
  }

   /**
   * The merged/consensus retention time end (latest rt) of this compound
   * @return rtEndSeconds
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_RT_END_SECONDS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public Double getRtEndSeconds() {
    return rtEndSeconds;
  }


  @JsonProperty(JSON_PROPERTY_RT_END_SECONDS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setRtEndSeconds(Double rtEndSeconds) {
    this.rtEndSeconds = rtEndSeconds;
  }

  public Compound neutralMass(Double neutralMass) {
    
    this.neutralMass = neutralMass;
    return this;
  }

   /**
   * Neutral mass of this compound. Ion masse minus the mass of the assigned adduct of each feature of  this compound should result in the same neutral mass
   * @return neutralMass
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_NEUTRAL_MASS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public Double getNeutralMass() {
    return neutralMass;
  }


  @JsonProperty(JSON_PROPERTY_NEUTRAL_MASS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setNeutralMass(Double neutralMass) {
    this.neutralMass = neutralMass;
  }

  public Compound features(List<AlignedFeature> features) {
    
    this.features = features;
    return this;
  }

  public Compound addFeaturesItem(AlignedFeature featuresItem) {
    if (this.features == null) {
      this.features = new ArrayList<>();
    }
    this.features.add(featuresItem);
    return this;
  }

   /**
   * List of aligned features (adducts) that belong to the same (this) compound
   * @return features
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_FEATURES)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public List<AlignedFeature> getFeatures() {
    return features;
  }


  @JsonProperty(JSON_PROPERTY_FEATURES)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setFeatures(List<AlignedFeature> features) {
    this.features = features;
  }

  public Compound consensusAnnotations(ConsensusAnnotationsCSI consensusAnnotations) {
    
    this.consensusAnnotations = consensusAnnotations;
    return this;
  }

   /**
   * Get consensusAnnotations
   * @return consensusAnnotations
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_CONSENSUS_ANNOTATIONS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public ConsensusAnnotationsCSI getConsensusAnnotations() {
    return consensusAnnotations;
  }


  @JsonProperty(JSON_PROPERTY_CONSENSUS_ANNOTATIONS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setConsensusAnnotations(ConsensusAnnotationsCSI consensusAnnotations) {
    this.consensusAnnotations = consensusAnnotations;
  }

  public Compound consensusAnnotationsDeNovo(ConsensusAnnotationsDeNovo consensusAnnotationsDeNovo) {
    
    this.consensusAnnotationsDeNovo = consensusAnnotationsDeNovo;
    return this;
  }

   /**
   * Get consensusAnnotationsDeNovo
   * @return consensusAnnotationsDeNovo
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_CONSENSUS_ANNOTATIONS_DE_NOVO)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public ConsensusAnnotationsDeNovo getConsensusAnnotationsDeNovo() {
    return consensusAnnotationsDeNovo;
  }


  @JsonProperty(JSON_PROPERTY_CONSENSUS_ANNOTATIONS_DE_NOVO)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setConsensusAnnotationsDeNovo(ConsensusAnnotationsDeNovo consensusAnnotationsDeNovo) {
    this.consensusAnnotationsDeNovo = consensusAnnotationsDeNovo;
  }

  public Compound customAnnotations(ConsensusAnnotationsCSI customAnnotations) {
    
    this.customAnnotations = customAnnotations;
    return this;
  }

   /**
   * Get customAnnotations
   * @return customAnnotations
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_CUSTOM_ANNOTATIONS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public ConsensusAnnotationsCSI getCustomAnnotations() {
    return customAnnotations;
  }


  @JsonProperty(JSON_PROPERTY_CUSTOM_ANNOTATIONS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setCustomAnnotations(ConsensusAnnotationsCSI customAnnotations) {
    this.customAnnotations = customAnnotations;
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    Compound compound = (Compound) o;
    return Objects.equals(this.compoundId, compound.compoundId) &&
        Objects.equals(this.name, compound.name) &&
        Objects.equals(this.rtStartSeconds, compound.rtStartSeconds) &&
        Objects.equals(this.rtEndSeconds, compound.rtEndSeconds) &&
        Objects.equals(this.neutralMass, compound.neutralMass) &&
        Objects.equals(this.features, compound.features) &&
        Objects.equals(this.consensusAnnotations, compound.consensusAnnotations) &&
        Objects.equals(this.consensusAnnotationsDeNovo, compound.consensusAnnotationsDeNovo) &&
        Objects.equals(this.customAnnotations, compound.customAnnotations);
  }

  @Override
  public int hashCode() {
    return Objects.hash(compoundId, name, rtStartSeconds, rtEndSeconds, neutralMass, features, consensusAnnotations, consensusAnnotationsDeNovo, customAnnotations);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class Compound {\n");
    sb.append("    compoundId: ").append(toIndentedString(compoundId)).append("\n");
    sb.append("    name: ").append(toIndentedString(name)).append("\n");
    sb.append("    rtStartSeconds: ").append(toIndentedString(rtStartSeconds)).append("\n");
    sb.append("    rtEndSeconds: ").append(toIndentedString(rtEndSeconds)).append("\n");
    sb.append("    neutralMass: ").append(toIndentedString(neutralMass)).append("\n");
    sb.append("    features: ").append(toIndentedString(features)).append("\n");
    sb.append("    consensusAnnotations: ").append(toIndentedString(consensusAnnotations)).append("\n");
    sb.append("    consensusAnnotationsDeNovo: ").append(toIndentedString(consensusAnnotationsDeNovo)).append("\n");
    sb.append("    customAnnotations: ").append(toIndentedString(customAnnotations)).append("\n");
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

