/*
 * SIRIUS Nightsky API
 * REST API that provides the full functionality of SIRIUS and its web services as background service. It is intended as entry-point for scripting languages and software integration SDKs.This API is exposed by SIRIUS 6
 *
 * The version of the OpenAPI document: 2.1
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
import de.unijena.bioinf.ms.nightsky.sdk.model.BasicSpectrum;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;
import com.fasterxml.jackson.annotation.JsonTypeName;

/**
 * IsotopePatternAnnotation
 */
@JsonPropertyOrder({
  IsotopePatternAnnotation.JSON_PROPERTY_ISOTOPE_PATTERN,
  IsotopePatternAnnotation.JSON_PROPERTY_SIMULATED_PATTERN
})
@jakarta.annotation.Generated(value = "org.openapitools.codegen.languages.JavaClientCodegen", comments = "Generator version: 7.6.0")
public class IsotopePatternAnnotation {
  public static final String JSON_PROPERTY_ISOTOPE_PATTERN = "isotopePattern";
  private BasicSpectrum isotopePattern;

  public static final String JSON_PROPERTY_SIMULATED_PATTERN = "simulatedPattern";
  private BasicSpectrum simulatedPattern;

  public IsotopePatternAnnotation() {
  }

  public IsotopePatternAnnotation isotopePattern(BasicSpectrum isotopePattern) {
    
    this.isotopePattern = isotopePattern;
    return this;
  }

   /**
   * Get isotopePattern
   * @return isotopePattern
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_ISOTOPE_PATTERN)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public BasicSpectrum getIsotopePattern() {
    return isotopePattern;
  }


  @JsonProperty(JSON_PROPERTY_ISOTOPE_PATTERN)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setIsotopePattern(BasicSpectrum isotopePattern) {
    this.isotopePattern = isotopePattern;
  }

  public IsotopePatternAnnotation simulatedPattern(BasicSpectrum simulatedPattern) {
    
    this.simulatedPattern = simulatedPattern;
    return this;
  }

   /**
   * Get simulatedPattern
   * @return simulatedPattern
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_SIMULATED_PATTERN)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public BasicSpectrum getSimulatedPattern() {
    return simulatedPattern;
  }


  @JsonProperty(JSON_PROPERTY_SIMULATED_PATTERN)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setSimulatedPattern(BasicSpectrum simulatedPattern) {
    this.simulatedPattern = simulatedPattern;
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    IsotopePatternAnnotation isotopePatternAnnotation = (IsotopePatternAnnotation) o;
    return Objects.equals(this.isotopePattern, isotopePatternAnnotation.isotopePattern) &&
        Objects.equals(this.simulatedPattern, isotopePatternAnnotation.simulatedPattern);
  }

  @Override
  public int hashCode() {
    return Objects.hash(isotopePattern, simulatedPattern);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class IsotopePatternAnnotation {\n");
    sb.append("    isotopePattern: ").append(toIndentedString(isotopePattern)).append("\n");
    sb.append("    simulatedPattern: ").append(toIndentedString(simulatedPattern)).append("\n");
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

