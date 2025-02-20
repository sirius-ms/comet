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
import com.fasterxml.jackson.annotation.JsonPropertyOrder;
import com.fasterxml.jackson.annotation.JsonTypeName;

/**
 * MsNovelist
 */
@JsonPropertyOrder({
  MsNovelist.JSON_PROPERTY_ENABLED,
  MsNovelist.JSON_PROPERTY_NUMBER_OF_CANDIDATE_TO_PREDICT
})
@jakarta.annotation.Generated(value = "org.openapitools.codegen.languages.JavaClientCodegen", comments = "Generator version: 7.6.0")
public class MsNovelist {
  public static final String JSON_PROPERTY_ENABLED = "enabled";
  private Boolean enabled;

  public static final String JSON_PROPERTY_NUMBER_OF_CANDIDATE_TO_PREDICT = "numberOfCandidateToPredict";
  private Integer numberOfCandidateToPredict;

  public MsNovelist() {
  }

  public MsNovelist enabled(Boolean enabled) {
    
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

  public MsNovelist numberOfCandidateToPredict(Integer numberOfCandidateToPredict) {
    
    this.numberOfCandidateToPredict = numberOfCandidateToPredict;
    return this;
  }

   /**
   * Number of structure candidates to be predicted by MsNovelist.  Max Value 128. Values &gt; 128 will be set to 128.  Actual number of returned candidate might be lower du to duplicates being created by MsNovelist.
   * @return numberOfCandidateToPredict
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_NUMBER_OF_CANDIDATE_TO_PREDICT)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public Integer getNumberOfCandidateToPredict() {
    return numberOfCandidateToPredict;
  }


  @JsonProperty(JSON_PROPERTY_NUMBER_OF_CANDIDATE_TO_PREDICT)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setNumberOfCandidateToPredict(Integer numberOfCandidateToPredict) {
    this.numberOfCandidateToPredict = numberOfCandidateToPredict;
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    MsNovelist msNovelist = (MsNovelist) o;
    return Objects.equals(this.enabled, msNovelist.enabled) &&
        Objects.equals(this.numberOfCandidateToPredict, msNovelist.numberOfCandidateToPredict);
  }

  @Override
  public int hashCode() {
    return Objects.hash(enabled, numberOfCandidateToPredict);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class MsNovelist {\n");
    sb.append("    enabled: ").append(toIndentedString(enabled)).append("\n");
    sb.append("    numberOfCandidateToPredict: ").append(toIndentedString(numberOfCandidateToPredict)).append("\n");
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

