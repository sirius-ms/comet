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
import de.unijena.bioinf.ms.nightsky.sdk.model.AnnotatedSpectrum;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;
import com.fasterxml.jackson.annotation.JsonTypeName;

/**
 * 
 */
@JsonPropertyOrder({
  IsotopePatternAnnotation.JSON_PROPERTY_ISOTOPE_PATTERN,
  IsotopePatternAnnotation.JSON_PROPERTY_SIMULATED_PATTERN,
  IsotopePatternAnnotation.JSON_PROPERTY_MS1,
  IsotopePatternAnnotation.JSON_PROPERTY_PEAKS_IN_MS1
})
@javax.annotation.Generated(value = "org.openapitools.codegen.languages.JavaClientCodegen")
public class IsotopePatternAnnotation {
  public static final String JSON_PROPERTY_ISOTOPE_PATTERN = "isotopePattern";
  private AnnotatedSpectrum isotopePattern;

  public static final String JSON_PROPERTY_SIMULATED_PATTERN = "simulatedPattern";
  private AnnotatedSpectrum simulatedPattern;

  public static final String JSON_PROPERTY_MS1 = "ms1";
  private AnnotatedSpectrum ms1;

  public static final String JSON_PROPERTY_PEAKS_IN_MS1 = "peaksInMs1";
  private List<Integer> peaksInMs1;

  public IsotopePatternAnnotation() {
  }

  public IsotopePatternAnnotation isotopePattern(AnnotatedSpectrum isotopePattern) {
    
    this.isotopePattern = isotopePattern;
    return this;
  }

   /**
   * Get isotopePattern
   * @return isotopePattern
  **/
  @javax.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_ISOTOPE_PATTERN)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public AnnotatedSpectrum getIsotopePattern() {
    return isotopePattern;
  }


  @JsonProperty(JSON_PROPERTY_ISOTOPE_PATTERN)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setIsotopePattern(AnnotatedSpectrum isotopePattern) {
    this.isotopePattern = isotopePattern;
  }


  public IsotopePatternAnnotation simulatedPattern(AnnotatedSpectrum simulatedPattern) {
    
    this.simulatedPattern = simulatedPattern;
    return this;
  }

   /**
   * Get simulatedPattern
   * @return simulatedPattern
  **/
  @javax.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_SIMULATED_PATTERN)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public AnnotatedSpectrum getSimulatedPattern() {
    return simulatedPattern;
  }


  @JsonProperty(JSON_PROPERTY_SIMULATED_PATTERN)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setSimulatedPattern(AnnotatedSpectrum simulatedPattern) {
    this.simulatedPattern = simulatedPattern;
  }


  public IsotopePatternAnnotation ms1(AnnotatedSpectrum ms1) {
    
    this.ms1 = ms1;
    return this;
  }

   /**
   * Get ms1
   * @return ms1
  **/
  @javax.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_MS1)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public AnnotatedSpectrum getMs1() {
    return ms1;
  }


  @JsonProperty(JSON_PROPERTY_MS1)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setMs1(AnnotatedSpectrum ms1) {
    this.ms1 = ms1;
  }


  public IsotopePatternAnnotation peaksInMs1(List<Integer> peaksInMs1) {
    
    this.peaksInMs1 = peaksInMs1;
    return this;
  }

  public IsotopePatternAnnotation addPeaksInMs1Item(Integer peaksInMs1Item) {
    if (this.peaksInMs1 == null) {
      this.peaksInMs1 = new ArrayList<>();
    }
    this.peaksInMs1.add(peaksInMs1Item);
    return this;
  }

   /**
   * Get peaksInMs1
   * @return peaksInMs1
  **/
  @javax.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_PEAKS_IN_MS1)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public List<Integer> getPeaksInMs1() {
    return peaksInMs1;
  }


  @JsonProperty(JSON_PROPERTY_PEAKS_IN_MS1)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setPeaksInMs1(List<Integer> peaksInMs1) {
    this.peaksInMs1 = peaksInMs1;
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
        Objects.equals(this.simulatedPattern, isotopePatternAnnotation.simulatedPattern) &&
        Objects.equals(this.ms1, isotopePatternAnnotation.ms1) &&
        Objects.equals(this.peaksInMs1, isotopePatternAnnotation.peaksInMs1);
  }

  @Override
  public int hashCode() {
    return Objects.hash(isotopePattern, simulatedPattern, ms1, peaksInMs1);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class IsotopePatternAnnotation {\n");
    sb.append("    isotopePattern: ").append(toIndentedString(isotopePattern)).append("\n");
    sb.append("    simulatedPattern: ").append(toIndentedString(simulatedPattern)).append("\n");
    sb.append("    ms1: ").append(toIndentedString(ms1)).append("\n");
    sb.append("    peaksInMs1: ").append(toIndentedString(peaksInMs1)).append("\n");
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

