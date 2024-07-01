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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;
import com.fasterxml.jackson.annotation.JsonTypeName;

/**
 * The MsData wraps all spectral input data belonging to a feature.  &lt;p&gt;  Each Feature has:  - One merged MS/MS spectrum (optional)  - One merged MS spectrum (optional)  - many MS/MS spectra  - many MS spectra  &lt;p&gt;  Each non-merged spectrum has an index which can be used to access the spectrum.  &lt;p&gt;  In the future we might add some additional information like chromatographic peak or something similar
 */
@JsonPropertyOrder({
  MsData.JSON_PROPERTY_MERGED_MS1,
  MsData.JSON_PROPERTY_MERGED_MS2,
  MsData.JSON_PROPERTY_MS1_SPECTRA,
  MsData.JSON_PROPERTY_MS2_SPECTRA
})
@jakarta.annotation.Generated(value = "org.openapitools.codegen.languages.JavaClientCodegen", comments = "Generator version: 7.6.0")
public class MsData {
  public static final String JSON_PROPERTY_MERGED_MS1 = "mergedMs1";
  private BasicSpectrum mergedMs1;

  public static final String JSON_PROPERTY_MERGED_MS2 = "mergedMs2";
  private BasicSpectrum mergedMs2;

  public static final String JSON_PROPERTY_MS1_SPECTRA = "ms1Spectra";
  private List<BasicSpectrum> ms1Spectra = new ArrayList<>();

  public static final String JSON_PROPERTY_MS2_SPECTRA = "ms2Spectra";
  private List<BasicSpectrum> ms2Spectra = new ArrayList<>();

  public MsData() {
  }

  public MsData mergedMs1(BasicSpectrum mergedMs1) {
    
    this.mergedMs1 = mergedMs1;
    return this;
  }

   /**
   * Get mergedMs1
   * @return mergedMs1
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_MERGED_MS1)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public BasicSpectrum getMergedMs1() {
    return mergedMs1;
  }


  @JsonProperty(JSON_PROPERTY_MERGED_MS1)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setMergedMs1(BasicSpectrum mergedMs1) {
    this.mergedMs1 = mergedMs1;
  }

  public MsData mergedMs2(BasicSpectrum mergedMs2) {
    
    this.mergedMs2 = mergedMs2;
    return this;
  }

   /**
   * Get mergedMs2
   * @return mergedMs2
  **/
  @jakarta.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_MERGED_MS2)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public BasicSpectrum getMergedMs2() {
    return mergedMs2;
  }


  @JsonProperty(JSON_PROPERTY_MERGED_MS2)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setMergedMs2(BasicSpectrum mergedMs2) {
    this.mergedMs2 = mergedMs2;
  }

  public MsData ms1Spectra(List<BasicSpectrum> ms1Spectra) {
    
    this.ms1Spectra = ms1Spectra;
    return this;
  }

  public MsData addMs1SpectraItem(BasicSpectrum ms1SpectraItem) {
    if (this.ms1Spectra == null) {
      this.ms1Spectra = new ArrayList<>();
    }
    this.ms1Spectra.add(ms1SpectraItem);
    return this;
  }

   /**
   * Get ms1Spectra
   * @return ms1Spectra
  **/
  @jakarta.annotation.Nonnull
  @JsonProperty(JSON_PROPERTY_MS1_SPECTRA)
  @JsonInclude(value = JsonInclude.Include.ALWAYS)

  public List<BasicSpectrum> getMs1Spectra() {
    return ms1Spectra;
  }


  @JsonProperty(JSON_PROPERTY_MS1_SPECTRA)
  @JsonInclude(value = JsonInclude.Include.ALWAYS)
  public void setMs1Spectra(List<BasicSpectrum> ms1Spectra) {
    this.ms1Spectra = ms1Spectra;
  }

  public MsData ms2Spectra(List<BasicSpectrum> ms2Spectra) {
    
    this.ms2Spectra = ms2Spectra;
    return this;
  }

  public MsData addMs2SpectraItem(BasicSpectrum ms2SpectraItem) {
    if (this.ms2Spectra == null) {
      this.ms2Spectra = new ArrayList<>();
    }
    this.ms2Spectra.add(ms2SpectraItem);
    return this;
  }

   /**
   * Get ms2Spectra
   * @return ms2Spectra
  **/
  @jakarta.annotation.Nonnull
  @JsonProperty(JSON_PROPERTY_MS2_SPECTRA)
  @JsonInclude(value = JsonInclude.Include.ALWAYS)

  public List<BasicSpectrum> getMs2Spectra() {
    return ms2Spectra;
  }


  @JsonProperty(JSON_PROPERTY_MS2_SPECTRA)
  @JsonInclude(value = JsonInclude.Include.ALWAYS)
  public void setMs2Spectra(List<BasicSpectrum> ms2Spectra) {
    this.ms2Spectra = ms2Spectra;
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    MsData msData = (MsData) o;
    return Objects.equals(this.mergedMs1, msData.mergedMs1) &&
        Objects.equals(this.mergedMs2, msData.mergedMs2) &&
        Objects.equals(this.ms1Spectra, msData.ms1Spectra) &&
        Objects.equals(this.ms2Spectra, msData.ms2Spectra);
  }

  @Override
  public int hashCode() {
    return Objects.hash(mergedMs1, mergedMs2, ms1Spectra, ms2Spectra);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class MsData {\n");
    sb.append("    mergedMs1: ").append(toIndentedString(mergedMs1)).append("\n");
    sb.append("    mergedMs2: ").append(toIndentedString(mergedMs2)).append("\n");
    sb.append("    ms1Spectra: ").append(toIndentedString(ms1Spectra)).append("\n");
    sb.append("    ms2Spectra: ").append(toIndentedString(ms2Spectra)).append("\n");
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

