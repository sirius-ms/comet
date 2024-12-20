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
import io.sirius.ms.sdk.model.AdductNetworkExperimental;
import io.sirius.ms.sdk.model.Axes;
import io.sirius.ms.sdk.model.TraceExperimental;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;
import com.fasterxml.jackson.annotation.JsonTypeName;

/**
 * EXPERIMENTAL: This schema is experimental and may be changed (or even removed) without notice until it is declared stable.
 */
@JsonPropertyOrder({
  TraceSetExperimental.JSON_PROPERTY_ADDUCT_NETWORK,
  TraceSetExperimental.JSON_PROPERTY_SAMPLE_ID,
  TraceSetExperimental.JSON_PROPERTY_SAMPLE_NAME,
  TraceSetExperimental.JSON_PROPERTY_AXES,
  TraceSetExperimental.JSON_PROPERTY_TRACES
})
@jakarta.annotation.Generated(value = "org.openapitools.codegen.languages.JavaClientCodegen", comments = "Generator version: 7.6.0")
public class TraceSetExperimental {
  public static final String JSON_PROPERTY_ADDUCT_NETWORK = "adductNetwork";
  private AdductNetworkExperimental adductNetwork;

  public static final String JSON_PROPERTY_SAMPLE_ID = "sampleId";
  private Long sampleId;

  public static final String JSON_PROPERTY_SAMPLE_NAME = "sampleName";
  private String sampleName;

  public static final String JSON_PROPERTY_AXES = "axes";
  private Axes axes;

  public static final String JSON_PROPERTY_TRACES = "traces";
  private List<TraceExperimental> traces = new ArrayList<>();

  public TraceSetExperimental() {
  }

  public TraceSetExperimental adductNetwork(AdductNetworkExperimental adductNetwork) {
    
    this.adductNetwork = adductNetwork;
    return this;
  }

   /**
   * Get adductNetwork
   * @return adductNetwork
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_ADDUCT_NETWORK)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public AdductNetworkExperimental getAdductNetwork() {
    return adductNetwork;
  }


  @JsonProperty(JSON_PROPERTY_ADDUCT_NETWORK)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setAdductNetwork(AdductNetworkExperimental adductNetwork) {
    this.adductNetwork = adductNetwork;
  }

  public TraceSetExperimental sampleId(Long sampleId) {
    
    this.sampleId = sampleId;
    return this;
  }

   /**
   * Get sampleId
   * @return sampleId
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_SAMPLE_ID)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public Long getSampleId() {
    return sampleId;
  }


  @JsonProperty(JSON_PROPERTY_SAMPLE_ID)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setSampleId(Long sampleId) {
    this.sampleId = sampleId;
  }

  public TraceSetExperimental sampleName(String sampleName) {
    
    this.sampleName = sampleName;
    return this;
  }

   /**
   * Get sampleName
   * @return sampleName
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_SAMPLE_NAME)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public String getSampleName() {
    return sampleName;
  }


  @JsonProperty(JSON_PROPERTY_SAMPLE_NAME)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setSampleName(String sampleName) {
    this.sampleName = sampleName;
  }

  public TraceSetExperimental axes(Axes axes) {
    
    this.axes = axes;
    return this;
  }

   /**
   * Get axes
   * @return axes
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_AXES)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public Axes getAxes() {
    return axes;
  }


  @JsonProperty(JSON_PROPERTY_AXES)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setAxes(Axes axes) {
    this.axes = axes;
  }

  public TraceSetExperimental traces(List<TraceExperimental> traces) {
    
    this.traces = traces;
    return this;
  }

  public TraceSetExperimental addTracesItem(TraceExperimental tracesItem) {
    if (this.traces == null) {
      this.traces = new ArrayList<>();
    }
    this.traces.add(tracesItem);
    return this;
  }

   /**
   * Get traces
   * @return traces
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_TRACES)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public List<TraceExperimental> getTraces() {
    return traces;
  }


  @JsonProperty(JSON_PROPERTY_TRACES)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setTraces(List<TraceExperimental> traces) {
    this.traces = traces;
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    TraceSetExperimental traceSetExperimental = (TraceSetExperimental) o;
    return Objects.equals(this.adductNetwork, traceSetExperimental.adductNetwork) &&
        Objects.equals(this.sampleId, traceSetExperimental.sampleId) &&
        Objects.equals(this.sampleName, traceSetExperimental.sampleName) &&
        Objects.equals(this.axes, traceSetExperimental.axes) &&
        Objects.equals(this.traces, traceSetExperimental.traces);
  }

  @Override
  public int hashCode() {
    return Objects.hash(adductNetwork, sampleId, sampleName, axes, traces);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class TraceSetExperimental {\n");
    sb.append("    adductNetwork: ").append(toIndentedString(adductNetwork)).append("\n");
    sb.append("    sampleId: ").append(toIndentedString(sampleId)).append("\n");
    sb.append("    sampleName: ").append(toIndentedString(sampleName)).append("\n");
    sb.append("    axes: ").append(toIndentedString(axes)).append("\n");
    sb.append("    traces: ").append(toIndentedString(traces)).append("\n");
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

