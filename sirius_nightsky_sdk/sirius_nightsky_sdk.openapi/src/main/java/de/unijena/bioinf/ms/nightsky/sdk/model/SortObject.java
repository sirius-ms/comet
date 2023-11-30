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
  SortObject.JSON_PROPERTY_EMPTY,
  SortObject.JSON_PROPERTY_SORTED,
  SortObject.JSON_PROPERTY_UNSORTED
})
@javax.annotation.Generated(value = "org.openapitools.codegen.languages.JavaClientCodegen")
public class SortObject {
  public static final String JSON_PROPERTY_EMPTY = "empty";
  private Boolean empty;

  public static final String JSON_PROPERTY_SORTED = "sorted";
  private Boolean sorted;

  public static final String JSON_PROPERTY_UNSORTED = "unsorted";
  private Boolean unsorted;

  public SortObject() {
  }

  public SortObject empty(Boolean empty) {
    
    this.empty = empty;
    return this;
  }

   /**
   * Get empty
   * @return empty
  **/
  @javax.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_EMPTY)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public Boolean isEmpty() {
    return empty;
  }


  @JsonProperty(JSON_PROPERTY_EMPTY)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setEmpty(Boolean empty) {
    this.empty = empty;
  }


  public SortObject sorted(Boolean sorted) {
    
    this.sorted = sorted;
    return this;
  }

   /**
   * Get sorted
   * @return sorted
  **/
  @javax.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_SORTED)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public Boolean isSorted() {
    return sorted;
  }


  @JsonProperty(JSON_PROPERTY_SORTED)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setSorted(Boolean sorted) {
    this.sorted = sorted;
  }


  public SortObject unsorted(Boolean unsorted) {
    
    this.unsorted = unsorted;
    return this;
  }

   /**
   * Get unsorted
   * @return unsorted
  **/
  @javax.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_UNSORTED)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public Boolean isUnsorted() {
    return unsorted;
  }


  @JsonProperty(JSON_PROPERTY_UNSORTED)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setUnsorted(Boolean unsorted) {
    this.unsorted = unsorted;
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    SortObject sortObject = (SortObject) o;
    return Objects.equals(this.empty, sortObject.empty) &&
        Objects.equals(this.sorted, sortObject.sorted) &&
        Objects.equals(this.unsorted, sortObject.unsorted);
  }

  @Override
  public int hashCode() {
    return Objects.hash(empty, sorted, unsorted);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class SortObject {\n");
    sb.append("    empty: ").append(toIndentedString(empty)).append("\n");
    sb.append("    sorted: ").append(toIndentedString(sorted)).append("\n");
    sb.append("    unsorted: ").append(toIndentedString(unsorted)).append("\n");
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

