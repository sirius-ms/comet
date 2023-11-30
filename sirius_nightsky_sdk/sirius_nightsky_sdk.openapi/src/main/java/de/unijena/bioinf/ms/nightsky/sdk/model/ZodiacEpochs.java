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
import com.fasterxml.jackson.annotation.JsonPropertyOrder;
import com.fasterxml.jackson.annotation.JsonTypeName;

/**
 * 
 */
@JsonPropertyOrder({
  ZodiacEpochs.JSON_PROPERTY_ITERATIONS,
  ZodiacEpochs.JSON_PROPERTY_BURN_IN_PERIOD,
  ZodiacEpochs.JSON_PROPERTY_NUMBER_OF_MARKOV_CHAINS
})
@javax.annotation.Generated(value = "org.openapitools.codegen.languages.JavaClientCodegen")
public class ZodiacEpochs {
  public static final String JSON_PROPERTY_ITERATIONS = "iterations";
  private Integer iterations;

  public static final String JSON_PROPERTY_BURN_IN_PERIOD = "burnInPeriod";
  private Integer burnInPeriod;

  public static final String JSON_PROPERTY_NUMBER_OF_MARKOV_CHAINS = "numberOfMarkovChains";
  private Integer numberOfMarkovChains;

  public ZodiacEpochs() {
  }

  public ZodiacEpochs iterations(Integer iterations) {
    
    this.iterations = iterations;
    return this;
  }

   /**
   * Get iterations
   * @return iterations
  **/
  @javax.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_ITERATIONS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public Integer getIterations() {
    return iterations;
  }


  @JsonProperty(JSON_PROPERTY_ITERATIONS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setIterations(Integer iterations) {
    this.iterations = iterations;
  }


  public ZodiacEpochs burnInPeriod(Integer burnInPeriod) {
    
    this.burnInPeriod = burnInPeriod;
    return this;
  }

   /**
   * Get burnInPeriod
   * @return burnInPeriod
  **/
  @javax.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_BURN_IN_PERIOD)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public Integer getBurnInPeriod() {
    return burnInPeriod;
  }


  @JsonProperty(JSON_PROPERTY_BURN_IN_PERIOD)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setBurnInPeriod(Integer burnInPeriod) {
    this.burnInPeriod = burnInPeriod;
  }


  public ZodiacEpochs numberOfMarkovChains(Integer numberOfMarkovChains) {
    
    this.numberOfMarkovChains = numberOfMarkovChains;
    return this;
  }

   /**
   * Get numberOfMarkovChains
   * @return numberOfMarkovChains
  **/
  @javax.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_NUMBER_OF_MARKOV_CHAINS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public Integer getNumberOfMarkovChains() {
    return numberOfMarkovChains;
  }


  @JsonProperty(JSON_PROPERTY_NUMBER_OF_MARKOV_CHAINS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setNumberOfMarkovChains(Integer numberOfMarkovChains) {
    this.numberOfMarkovChains = numberOfMarkovChains;
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    ZodiacEpochs zodiacEpochs = (ZodiacEpochs) o;
    return Objects.equals(this.iterations, zodiacEpochs.iterations) &&
        Objects.equals(this.burnInPeriod, zodiacEpochs.burnInPeriod) &&
        Objects.equals(this.numberOfMarkovChains, zodiacEpochs.numberOfMarkovChains);
  }

  @Override
  public int hashCode() {
    return Objects.hash(iterations, burnInPeriod, numberOfMarkovChains);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class ZodiacEpochs {\n");
    sb.append("    iterations: ").append(toIndentedString(iterations)).append("\n");
    sb.append("    burnInPeriod: ").append(toIndentedString(burnInPeriod)).append("\n");
    sb.append("    numberOfMarkovChains: ").append(toIndentedString(numberOfMarkovChains)).append("\n");
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

