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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;
import com.fasterxml.jackson.annotation.JsonTypeName;

/**
 * WorkerInfo
 */
@JsonPropertyOrder({
  WorkerInfo.JSON_PROPERTY_ID,
  WorkerInfo.JSON_PROPERTY_TYPE,
  WorkerInfo.JSON_PROPERTY_SUPPORTED_PREDICTORS,
  WorkerInfo.JSON_PROPERTY_VERSION,
  WorkerInfo.JSON_PROPERTY_HOST,
  WorkerInfo.JSON_PROPERTY_PREFIX,
  WorkerInfo.JSON_PROPERTY_STATE,
  WorkerInfo.JSON_PROPERTY_ALIVE,
  WorkerInfo.JSON_PROPERTY_SERVER_TIME
})
@jakarta.annotation.Generated(value = "org.openapitools.codegen.languages.JavaClientCodegen", comments = "Generator version: 7.6.0")
public class WorkerInfo {
  public static final String JSON_PROPERTY_ID = "id";
  private Long id;

  /**
   * Gets or Sets type
   */
  public enum TypeEnum {
    FORMULA_ID("FORMULA_ID"),
    
    FINGER_ID("FINGER_ID"),
    
    IOKR("IOKR"),
    
    CANOPUS("CANOPUS"),
    
    COVTREE("COVTREE");

    private String value;

    TypeEnum(String value) {
      this.value = value;
    }

    @JsonValue
    public String getValue() {
      return value;
    }

    @Override
    public String toString() {
      return String.valueOf(value);
    }

    @JsonCreator
    public static TypeEnum fromValue(String value) {
      for (TypeEnum b : TypeEnum.values()) {
        if (b.value.equals(value)) {
          return b;
        }
      }
      throw new IllegalArgumentException("Unexpected value '" + value + "'");
    }
  }

  public static final String JSON_PROPERTY_TYPE = "type";
  private TypeEnum type;

  /**
   * Gets or Sets supportedPredictors
   */
  public enum SupportedPredictorsEnum {
    POSITIVE("CSI_FINGERID_POSITIVE"),
    
    NEGATIVE("CSI_FINGERID_NEGATIVE");

    private String value;

    SupportedPredictorsEnum(String value) {
      this.value = value;
    }

    @JsonValue
    public String getValue() {
      return value;
    }

    @Override
    public String toString() {
      return String.valueOf(value);
    }

    @JsonCreator
    public static SupportedPredictorsEnum fromValue(String value) {
      for (SupportedPredictorsEnum b : SupportedPredictorsEnum.values()) {
        if (b.value.equals(value)) {
          return b;
        }
      }
      throw new IllegalArgumentException("Unexpected value '" + value + "'");
    }
  }

  public static final String JSON_PROPERTY_SUPPORTED_PREDICTORS = "supportedPredictors";
  private List<SupportedPredictorsEnum> supportedPredictors = new ArrayList<>();

  public static final String JSON_PROPERTY_VERSION = "version";
  private String version;

  public static final String JSON_PROPERTY_HOST = "host";
  private String host;

  public static final String JSON_PROPERTY_PREFIX = "prefix";
  private String prefix;

  public static final String JSON_PROPERTY_STATE = "state";
  private Integer state;

  public static final String JSON_PROPERTY_ALIVE = "alive";
  private Long alive;

  public static final String JSON_PROPERTY_SERVER_TIME = "serverTime";
  private Long serverTime;

  public WorkerInfo() {
  }

  public WorkerInfo id(Long id) {
    
    this.id = id;
    return this;
  }

   /**
   * Get id
   * @return id
  **/
  @jakarta.annotation.Nonnull
  @JsonProperty(JSON_PROPERTY_ID)
  @JsonInclude(value = JsonInclude.Include.ALWAYS)

  public Long getId() {
    return id;
  }


  @JsonProperty(JSON_PROPERTY_ID)
  @JsonInclude(value = JsonInclude.Include.ALWAYS)
  public void setId(Long id) {
    this.id = id;
  }

  public WorkerInfo type(TypeEnum type) {
    
    this.type = type;
    return this;
  }

   /**
   * Get type
   * @return type
  **/
  @jakarta.annotation.Nonnull
  @JsonProperty(JSON_PROPERTY_TYPE)
  @JsonInclude(value = JsonInclude.Include.ALWAYS)

  public TypeEnum getType() {
    return type;
  }


  @JsonProperty(JSON_PROPERTY_TYPE)
  @JsonInclude(value = JsonInclude.Include.ALWAYS)
  public void setType(TypeEnum type) {
    this.type = type;
  }

  public WorkerInfo supportedPredictors(List<SupportedPredictorsEnum> supportedPredictors) {
    
    this.supportedPredictors = supportedPredictors;
    return this;
  }

  public WorkerInfo addSupportedPredictorsItem(SupportedPredictorsEnum supportedPredictorsItem) {
    if (this.supportedPredictors == null) {
      this.supportedPredictors = new ArrayList<>();
    }
    this.supportedPredictors.add(supportedPredictorsItem);
    return this;
  }

   /**
   * Get supportedPredictors
   * @return supportedPredictors
  **/
  @jakarta.annotation.Nonnull
  @JsonProperty(JSON_PROPERTY_SUPPORTED_PREDICTORS)
  @JsonInclude(value = JsonInclude.Include.ALWAYS)

  public List<SupportedPredictorsEnum> getSupportedPredictors() {
    return supportedPredictors;
  }


  @JsonProperty(JSON_PROPERTY_SUPPORTED_PREDICTORS)
  @JsonInclude(value = JsonInclude.Include.ALWAYS)
  public void setSupportedPredictors(List<SupportedPredictorsEnum> supportedPredictors) {
    this.supportedPredictors = supportedPredictors;
  }

  public WorkerInfo version(String version) {
    
    this.version = version;
    return this;
  }

   /**
   * Get version
   * @return version
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_VERSION)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public String getVersion() {
    return version;
  }


  @JsonProperty(JSON_PROPERTY_VERSION)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setVersion(String version) {
    this.version = version;
  }

  public WorkerInfo host(String host) {
    
    this.host = host;
    return this;
  }

   /**
   * Get host
   * @return host
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_HOST)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public String getHost() {
    return host;
  }


  @JsonProperty(JSON_PROPERTY_HOST)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setHost(String host) {
    this.host = host;
  }

  public WorkerInfo prefix(String prefix) {
    
    this.prefix = prefix;
    return this;
  }

   /**
   * Get prefix
   * @return prefix
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_PREFIX)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public String getPrefix() {
    return prefix;
  }


  @JsonProperty(JSON_PROPERTY_PREFIX)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setPrefix(String prefix) {
    this.prefix = prefix;
  }

  public WorkerInfo state(Integer state) {
    
    this.state = state;
    return this;
  }

   /**
   * Get state
   * @return state
  **/
  @jakarta.annotation.Nonnull
  @JsonProperty(JSON_PROPERTY_STATE)
  @JsonInclude(value = JsonInclude.Include.ALWAYS)

  public Integer getState() {
    return state;
  }


  @JsonProperty(JSON_PROPERTY_STATE)
  @JsonInclude(value = JsonInclude.Include.ALWAYS)
  public void setState(Integer state) {
    this.state = state;
  }

  public WorkerInfo alive(Long alive) {
    
    this.alive = alive;
    return this;
  }

   /**
   * Get alive
   * @return alive
  **/
  @jakarta.annotation.Nonnull
  @JsonProperty(JSON_PROPERTY_ALIVE)
  @JsonInclude(value = JsonInclude.Include.ALWAYS)

  public Long getAlive() {
    return alive;
  }


  @JsonProperty(JSON_PROPERTY_ALIVE)
  @JsonInclude(value = JsonInclude.Include.ALWAYS)
  public void setAlive(Long alive) {
    this.alive = alive;
  }

  public WorkerInfo serverTime(Long serverTime) {
    
    this.serverTime = serverTime;
    return this;
  }

   /**
   * Get serverTime
   * @return serverTime
  **/
  @jakarta.annotation.Nonnull
  @JsonProperty(JSON_PROPERTY_SERVER_TIME)
  @JsonInclude(value = JsonInclude.Include.ALWAYS)

  public Long getServerTime() {
    return serverTime;
  }


  @JsonProperty(JSON_PROPERTY_SERVER_TIME)
  @JsonInclude(value = JsonInclude.Include.ALWAYS)
  public void setServerTime(Long serverTime) {
    this.serverTime = serverTime;
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    WorkerInfo workerInfo = (WorkerInfo) o;
    return Objects.equals(this.id, workerInfo.id) &&
        Objects.equals(this.type, workerInfo.type) &&
        Objects.equals(this.supportedPredictors, workerInfo.supportedPredictors) &&
        Objects.equals(this.version, workerInfo.version) &&
        Objects.equals(this.host, workerInfo.host) &&
        Objects.equals(this.prefix, workerInfo.prefix) &&
        Objects.equals(this.state, workerInfo.state) &&
        Objects.equals(this.alive, workerInfo.alive) &&
        Objects.equals(this.serverTime, workerInfo.serverTime);
  }

  @Override
  public int hashCode() {
    return Objects.hash(id, type, supportedPredictors, version, host, prefix, state, alive, serverTime);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class WorkerInfo {\n");
    sb.append("    id: ").append(toIndentedString(id)).append("\n");
    sb.append("    type: ").append(toIndentedString(type)).append("\n");
    sb.append("    supportedPredictors: ").append(toIndentedString(supportedPredictors)).append("\n");
    sb.append("    version: ").append(toIndentedString(version)).append("\n");
    sb.append("    host: ").append(toIndentedString(host)).append("\n");
    sb.append("    prefix: ").append(toIndentedString(prefix)).append("\n");
    sb.append("    state: ").append(toIndentedString(state)).append("\n");
    sb.append("    alive: ").append(toIndentedString(alive)).append("\n");
    sb.append("    serverTime: ").append(toIndentedString(serverTime)).append("\n");
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

