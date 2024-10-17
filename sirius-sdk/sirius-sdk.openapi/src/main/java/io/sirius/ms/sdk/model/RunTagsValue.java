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


package io.sirius.ms.sdk.model;

import java.util.Objects;
import java.util.Arrays;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonTypeName;
import com.fasterxml.jackson.annotation.JsonValue;
import io.sirius.ms.sdk.model.BoolTag;
import io.sirius.ms.sdk.model.DoubleTag;
import io.sirius.ms.sdk.model.IntTag;
import io.sirius.ms.sdk.model.StringTag;
import io.sirius.ms.sdk.model.Tag;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;
import com.fasterxml.jackson.annotation.JsonTypeName;

/**
 * RunTagsValue
 */
@JsonPropertyOrder({
  RunTagsValue.JSON_PROPERTY_CATEGORY_NAME,
  RunTagsValue.JSON_PROPERTY_VALUE
})
@JsonTypeName("Run_tags_value")
@jakarta.annotation.Generated(value = "org.openapitools.codegen.languages.JavaClientCodegen", comments = "Generator version: 7.6.0")
public class RunTagsValue {
  public static final String JSON_PROPERTY_CATEGORY_NAME = "categoryName";
  private String categoryName;

  public static final String JSON_PROPERTY_VALUE = "value";
  private String value;

  public RunTagsValue() {
  }

  public RunTagsValue categoryName(String categoryName) {
    
    this.categoryName = categoryName;
    return this;
  }

   /**
   * Name of the tag category
   * @return categoryName
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_CATEGORY_NAME)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public String getCategoryName() {
    return categoryName;
  }


  @JsonProperty(JSON_PROPERTY_CATEGORY_NAME)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setCategoryName(String categoryName) {
    this.categoryName = categoryName;
  }

  public RunTagsValue value(String value) {
    
    this.value = value;
    return this;
  }

   /**
   * Tag value
   * @return value
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_VALUE)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public String getValue() {
    return value;
  }


  @JsonProperty(JSON_PROPERTY_VALUE)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setValue(String value) {
    this.value = value;
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    RunTagsValue runTagsValue = (RunTagsValue) o;
    return Objects.equals(this.categoryName, runTagsValue.categoryName) &&
        Objects.equals(this.value, runTagsValue.value);
  }

  @Override
  public int hashCode() {
    return Objects.hash(categoryName, value);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class RunTagsValue {\n");
    sb.append("    categoryName: ").append(toIndentedString(categoryName)).append("\n");
    sb.append("    value: ").append(toIndentedString(value)).append("\n");
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

