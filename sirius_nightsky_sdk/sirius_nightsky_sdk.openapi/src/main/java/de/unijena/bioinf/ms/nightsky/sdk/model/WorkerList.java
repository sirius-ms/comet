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

import java.net.URLEncoder;
import java.nio.charset.StandardCharsets;
import java.util.StringJoiner;
import java.util.Objects;
import java.util.Map;
import java.util.HashMap;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonTypeName;
import com.fasterxml.jackson.annotation.JsonValue;
import de.unijena.bioinf.ms.nightsky.sdk.model.WorkerInfo;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * 
 */
@JsonPropertyOrder({
  WorkerList.JSON_PROPERTY_PENDING_JOBS,
  WorkerList.JSON_PROPERTY_WORKER_LIST
})
@javax.annotation.Generated(value = "org.openapitools.codegen.languages.JavaClientCodegen")
public class WorkerList {
  public static final String JSON_PROPERTY_PENDING_JOBS = "pendingJobs";
  private Integer pendingJobs;

  public static final String JSON_PROPERTY_WORKER_LIST = "workerList";
  private List<WorkerInfo> workerList = new ArrayList<>();

  public WorkerList() { 
  }

  public WorkerList pendingJobs(Integer pendingJobs) {
    this.pendingJobs = pendingJobs;
    return this;
  }

   /**
   * Get pendingJobs
   * @return pendingJobs
  **/
  @javax.annotation.Nonnull
  @JsonProperty(JSON_PROPERTY_PENDING_JOBS)
  @JsonInclude(value = JsonInclude.Include.ALWAYS)

  public Integer getPendingJobs() {
    return pendingJobs;
  }


  @JsonProperty(JSON_PROPERTY_PENDING_JOBS)
  @JsonInclude(value = JsonInclude.Include.ALWAYS)
  public void setPendingJobs(Integer pendingJobs) {
    this.pendingJobs = pendingJobs;
  }


  public WorkerList workerList(List<WorkerInfo> workerList) {
    this.workerList = workerList;
    return this;
  }

  public WorkerList addWorkerListItem(WorkerInfo workerListItem) {
    if (this.workerList == null) {
      this.workerList = new ArrayList<>();
    }
    this.workerList.add(workerListItem);
    return this;
  }

   /**
   * Get workerList
   * @return workerList
  **/
  @javax.annotation.Nonnull
  @JsonProperty(JSON_PROPERTY_WORKER_LIST)
  @JsonInclude(value = JsonInclude.Include.ALWAYS)

  public List<WorkerInfo> getWorkerList() {
    return workerList;
  }


  @JsonProperty(JSON_PROPERTY_WORKER_LIST)
  @JsonInclude(value = JsonInclude.Include.ALWAYS)
  public void setWorkerList(List<WorkerInfo> workerList) {
    this.workerList = workerList;
  }


  /**
   * Return true if this WorkerList object is equal to o.
   */
  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    WorkerList workerList = (WorkerList) o;
    return Objects.equals(this.pendingJobs, workerList.pendingJobs) &&
        Objects.equals(this.workerList, workerList.workerList);
  }

  @Override
  public int hashCode() {
    return Objects.hash(pendingJobs, workerList);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class WorkerList {\n");
    sb.append("    pendingJobs: ").append(toIndentedString(pendingJobs)).append("\n");
    sb.append("    workerList: ").append(toIndentedString(workerList)).append("\n");
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

  /**
   * Convert the instance into URL query string.
   *
   * @return URL query string
   */
  public String toUrlQueryString() {
    return toUrlQueryString(null);
  }

  /**
   * Convert the instance into URL query string.
   *
   * @param prefix prefix of the query string
   * @return URL query string
   */
  public String toUrlQueryString(String prefix) {
    String suffix = "";
    String containerSuffix = "";
    String containerPrefix = "";
    if (prefix == null) {
      // style=form, explode=true, e.g. /pet?name=cat&type=manx
      prefix = "";
    } else {
      // deepObject style e.g. /pet?id[name]=cat&id[type]=manx
      prefix = prefix + "[";
      suffix = "]";
      containerSuffix = "]";
      containerPrefix = "[";
    }

    StringJoiner joiner = new StringJoiner("&");

    // add `pendingJobs` to the URL query string
    if (getPendingJobs() != null) {
      joiner.add(String.format("%spendingJobs%s=%s", prefix, suffix, URLEncoder.encode(String.valueOf(getPendingJobs()), StandardCharsets.UTF_8).replaceAll("\\+", "%20")));
    }

    // add `workerList` to the URL query string
    if (getWorkerList() != null) {
      for (int i = 0; i < getWorkerList().size(); i++) {
        if (getWorkerList().get(i) != null) {
          joiner.add(getWorkerList().get(i).toUrlQueryString(String.format("%sworkerList%s%s", prefix, suffix,
          "".equals(suffix) ? "" : String.format("%s%d%s", containerPrefix, i, containerSuffix))));
        }
      }
    }

    return joiner.toString();
  }
}

