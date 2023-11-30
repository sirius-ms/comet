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
import de.unijena.bioinf.ms.nightsky.sdk.model.DBLink;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;
import com.fasterxml.jackson.annotation.JsonTypeName;

/**
 * 
 */
@JsonPropertyOrder({
  StructureCandidate.JSON_PROPERTY_INCHI_KEY,
  StructureCandidate.JSON_PROPERTY_SMILES,
  StructureCandidate.JSON_PROPERTY_STRUCTURE_NAME,
  StructureCandidate.JSON_PROPERTY_XLOG_P,
  StructureCandidate.JSON_PROPERTY_DB_LINKS,
  StructureCandidate.JSON_PROPERTY_REF_SPECTRA_LINKS
})
@javax.annotation.Generated(value = "org.openapitools.codegen.languages.JavaClientCodegen")
public class StructureCandidate {
  public static final String JSON_PROPERTY_INCHI_KEY = "inchiKey";
  private String inchiKey;

  public static final String JSON_PROPERTY_SMILES = "smiles";
  private String smiles;

  public static final String JSON_PROPERTY_STRUCTURE_NAME = "structureName";
  private String structureName;

  public static final String JSON_PROPERTY_XLOG_P = "xlogP";
  private Double xlogP;

  public static final String JSON_PROPERTY_DB_LINKS = "dbLinks";
  private List<DBLink> dbLinks;

  public static final String JSON_PROPERTY_REF_SPECTRA_LINKS = "refSpectraLinks";
  private List<DBLink> refSpectraLinks;

  public StructureCandidate() {
  }

  public StructureCandidate inchiKey(String inchiKey) {
    
    this.inchiKey = inchiKey;
    return this;
  }

   /**
   * Get inchiKey
   * @return inchiKey
  **/
  @javax.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_INCHI_KEY)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public String getInchiKey() {
    return inchiKey;
  }


  @JsonProperty(JSON_PROPERTY_INCHI_KEY)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setInchiKey(String inchiKey) {
    this.inchiKey = inchiKey;
  }


  public StructureCandidate smiles(String smiles) {
    
    this.smiles = smiles;
    return this;
  }

   /**
   * Get smiles
   * @return smiles
  **/
  @javax.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_SMILES)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public String getSmiles() {
    return smiles;
  }


  @JsonProperty(JSON_PROPERTY_SMILES)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setSmiles(String smiles) {
    this.smiles = smiles;
  }


  public StructureCandidate structureName(String structureName) {
    
    this.structureName = structureName;
    return this;
  }

   /**
   * Get structureName
   * @return structureName
  **/
  @javax.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_STRUCTURE_NAME)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public String getStructureName() {
    return structureName;
  }


  @JsonProperty(JSON_PROPERTY_STRUCTURE_NAME)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setStructureName(String structureName) {
    this.structureName = structureName;
  }


  public StructureCandidate xlogP(Double xlogP) {
    
    this.xlogP = xlogP;
    return this;
  }

   /**
   * Get xlogP
   * @return xlogP
  **/
  @javax.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_XLOG_P)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public Double getXlogP() {
    return xlogP;
  }


  @JsonProperty(JSON_PROPERTY_XLOG_P)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setXlogP(Double xlogP) {
    this.xlogP = xlogP;
  }


  public StructureCandidate dbLinks(List<DBLink> dbLinks) {
    
    this.dbLinks = dbLinks;
    return this;
  }

  public StructureCandidate addDbLinksItem(DBLink dbLinksItem) {
    if (this.dbLinks == null) {
      this.dbLinks = new ArrayList<>();
    }
    this.dbLinks.add(dbLinksItem);
    return this;
  }

   /**
   * List of structure database links belonging to this structure candidate  OPTIONAL: needs to be added by parameter
   * @return dbLinks
  **/
  @javax.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_DB_LINKS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public List<DBLink> getDbLinks() {
    return dbLinks;
  }


  @JsonProperty(JSON_PROPERTY_DB_LINKS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setDbLinks(List<DBLink> dbLinks) {
    this.dbLinks = dbLinks;
  }


  public StructureCandidate refSpectraLinks(List<DBLink> refSpectraLinks) {
    
    this.refSpectraLinks = refSpectraLinks;
    return this;
  }

  public StructureCandidate addRefSpectraLinksItem(DBLink refSpectraLinksItem) {
    if (this.refSpectraLinks == null) {
      this.refSpectraLinks = new ArrayList<>();
    }
    this.refSpectraLinks.add(refSpectraLinksItem);
    return this;
  }

   /**
   * List of spectral library links belonging to this structure candidate  OPTIONAL: needs to be added by parameter
   * @return refSpectraLinks
  **/
  @javax.annotation.Nullable
  @JsonProperty(JSON_PROPERTY_REF_SPECTRA_LINKS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)

  public List<DBLink> getRefSpectraLinks() {
    return refSpectraLinks;
  }


  @JsonProperty(JSON_PROPERTY_REF_SPECTRA_LINKS)
  @JsonInclude(value = JsonInclude.Include.USE_DEFAULTS)
  public void setRefSpectraLinks(List<DBLink> refSpectraLinks) {
    this.refSpectraLinks = refSpectraLinks;
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    StructureCandidate structureCandidate = (StructureCandidate) o;
    return Objects.equals(this.inchiKey, structureCandidate.inchiKey) &&
        Objects.equals(this.smiles, structureCandidate.smiles) &&
        Objects.equals(this.structureName, structureCandidate.structureName) &&
        Objects.equals(this.xlogP, structureCandidate.xlogP) &&
        Objects.equals(this.dbLinks, structureCandidate.dbLinks) &&
        Objects.equals(this.refSpectraLinks, structureCandidate.refSpectraLinks);
  }

  @Override
  public int hashCode() {
    return Objects.hash(inchiKey, smiles, structureName, xlogP, dbLinks, refSpectraLinks);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class StructureCandidate {\n");
    sb.append("    inchiKey: ").append(toIndentedString(inchiKey)).append("\n");
    sb.append("    smiles: ").append(toIndentedString(smiles)).append("\n");
    sb.append("    structureName: ").append(toIndentedString(structureName)).append("\n");
    sb.append("    xlogP: ").append(toIndentedString(xlogP)).append("\n");
    sb.append("    dbLinks: ").append(toIndentedString(dbLinks)).append("\n");
    sb.append("    refSpectraLinks: ").append(toIndentedString(refSpectraLinks)).append("\n");
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

