package de.unijena.bioinf.ms.nightsky.sdk.api;

import de.unijena.bioinf.ms.nightsky.sdk.client.ApiClient;

import de.unijena.bioinf.ms.nightsky.sdk.model.GuiParameters;

import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.stream.Collectors;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.util.LinkedMultiValueMap;
import org.springframework.util.MultiValueMap;
import org.springframework.core.ParameterizedTypeReference;
import org.springframework.web.reactive.function.client.WebClient.ResponseSpec;
import org.springframework.web.reactive.function.client.WebClientResponseException;
import org.springframework.core.io.FileSystemResource;
import org.springframework.http.HttpHeaders;
import org.springframework.http.HttpMethod;
import org.springframework.http.HttpStatus;
import org.springframework.http.MediaType;
import org.springframework.http.ResponseEntity;
import reactor.core.publisher.Mono;
import reactor.core.publisher.Flux;

@javax.annotation.Generated(value = "org.openapitools.codegen.languages.JavaClientCodegen")
public class ExperimentalGuiApi {
    private ApiClient apiClient;

    public ExperimentalGuiApi() {
        this(new ApiClient());
    }

    @Autowired
    public ExperimentalGuiApi(ApiClient apiClient) {
        this.apiClient = apiClient;
    }

    public ApiClient getApiClient() {
        return apiClient;
    }

    public void setApiClient(ApiClient apiClient) {
        this.apiClient = apiClient;
    }

    /**
     * Apply given changes to the running GUI instance.
     * Apply given changes to the running GUI instance.
     * <p><b>200</b> - OK
     * @param projectId of project-space the GUI instance is connected to.
     * @param guiParameters parameters that should be applied.
     * @throws WebClientResponseException if an error occurs while attempting to invoke the API
     */
    private ResponseSpec applyToGuiRequestCreation(String projectId, GuiParameters guiParameters) throws WebClientResponseException {
        Object postBody = guiParameters;
        // verify the required parameter 'projectId' is set
        if (projectId == null) {
            throw new WebClientResponseException("Missing the required parameter 'projectId' when calling applyToGui", HttpStatus.BAD_REQUEST.value(), HttpStatus.BAD_REQUEST.getReasonPhrase(), null, null, null);
        }
        // verify the required parameter 'guiParameters' is set
        if (guiParameters == null) {
            throw new WebClientResponseException("Missing the required parameter 'guiParameters' when calling applyToGui", HttpStatus.BAD_REQUEST.value(), HttpStatus.BAD_REQUEST.getReasonPhrase(), null, null, null);
        }
        // create path and map variables
        final Map<String, Object> pathParams = new HashMap<String, Object>();

        pathParams.put("projectId", projectId);

        final MultiValueMap<String, String> queryParams = new LinkedMultiValueMap<String, String>();
        final HttpHeaders headerParams = new HttpHeaders();
        final MultiValueMap<String, String> cookieParams = new LinkedMultiValueMap<String, String>();
        final MultiValueMap<String, Object> formParams = new LinkedMultiValueMap<String, Object>();

        final String[] localVarAccepts = { };
        final List<MediaType> localVarAccept = apiClient.selectHeaderAccept(localVarAccepts);
        final String[] localVarContentTypes = { 
            "application/json"
        };
        final MediaType localVarContentType = apiClient.selectHeaderContentType(localVarContentTypes);

        String[] localVarAuthNames = new String[] {  };

        ParameterizedTypeReference<Void> localVarReturnType = new ParameterizedTypeReference<Void>() {};
        return apiClient.invokeAPI("/api/projects/{projectId}/gui", HttpMethod.PATCH, pathParams, queryParams, postBody, headerParams, cookieParams, formParams, localVarAccept, localVarContentType, localVarAuthNames, localVarReturnType);
    }

    /**
     * Apply given changes to the running GUI instance.
     * Apply given changes to the running GUI instance.
     * <p><b>200</b> - OK
     * @param projectId of project-space the GUI instance is connected to.
     * @param guiParameters parameters that should be applied.
     * @throws WebClientResponseException if an error occurs while attempting to invoke the API
     */
    public void applyToGui(String projectId, GuiParameters guiParameters) throws WebClientResponseException {
        ParameterizedTypeReference<Void> localVarReturnType = new ParameterizedTypeReference<Void>() {};
        applyToGuiRequestCreation(projectId, guiParameters).bodyToMono(localVarReturnType).block();
    }

    /**
     * Apply given changes to the running GUI instance.
     * Apply given changes to the running GUI instance.
     * <p><b>200</b> - OK
     * @param projectId of project-space the GUI instance is connected to.
     * @param guiParameters parameters that should be applied.
     * @throws WebClientResponseException if an error occurs while attempting to invoke the API
     */
    public ResponseEntity<Void> applyToGuiWithHttpInfo(String projectId, GuiParameters guiParameters) throws WebClientResponseException {
        ParameterizedTypeReference<Void> localVarReturnType = new ParameterizedTypeReference<Void>() {};
        return applyToGuiRequestCreation(projectId, guiParameters).toEntity(localVarReturnType).block();
    }

    /**
     * Apply given changes to the running GUI instance.
     * Apply given changes to the running GUI instance.
     * <p><b>200</b> - OK
     * @param projectId of project-space the GUI instance is connected to.
     * @param guiParameters parameters that should be applied.
     * @return ResponseSpec
     * @throws WebClientResponseException if an error occurs while attempting to invoke the API
     */
    public ResponseSpec applyToGuiWithResponseSpec(String projectId, GuiParameters guiParameters) throws WebClientResponseException {
        return applyToGuiRequestCreation(projectId, guiParameters);
    }
    /**
     * Close GUI instance of given project-space if available.
     * Close GUI instance of given project-space if available.
     * <p><b>200</b> - OK
     * @param projectId if project-space the GUI instance is connected to.
     * @throws WebClientResponseException if an error occurs while attempting to invoke the API
     */
    private ResponseSpec closeGuiRequestCreation(String projectId) throws WebClientResponseException {
        Object postBody = null;
        // verify the required parameter 'projectId' is set
        if (projectId == null) {
            throw new WebClientResponseException("Missing the required parameter 'projectId' when calling closeGui", HttpStatus.BAD_REQUEST.value(), HttpStatus.BAD_REQUEST.getReasonPhrase(), null, null, null);
        }
        // create path and map variables
        final Map<String, Object> pathParams = new HashMap<String, Object>();

        pathParams.put("projectId", projectId);

        final MultiValueMap<String, String> queryParams = new LinkedMultiValueMap<String, String>();
        final HttpHeaders headerParams = new HttpHeaders();
        final MultiValueMap<String, String> cookieParams = new LinkedMultiValueMap<String, String>();
        final MultiValueMap<String, Object> formParams = new LinkedMultiValueMap<String, Object>();

        final String[] localVarAccepts = { };
        final List<MediaType> localVarAccept = apiClient.selectHeaderAccept(localVarAccepts);
        final String[] localVarContentTypes = { };
        final MediaType localVarContentType = apiClient.selectHeaderContentType(localVarContentTypes);

        String[] localVarAuthNames = new String[] {  };

        ParameterizedTypeReference<Void> localVarReturnType = new ParameterizedTypeReference<Void>() {};
        return apiClient.invokeAPI("/api/projects/{projectId}/gui", HttpMethod.DELETE, pathParams, queryParams, postBody, headerParams, cookieParams, formParams, localVarAccept, localVarContentType, localVarAuthNames, localVarReturnType);
    }

    /**
     * Close GUI instance of given project-space if available.
     * Close GUI instance of given project-space if available.
     * <p><b>200</b> - OK
     * @param projectId if project-space the GUI instance is connected to.
     * @throws WebClientResponseException if an error occurs while attempting to invoke the API
     */
    public void closeGui(String projectId) throws WebClientResponseException {
        ParameterizedTypeReference<Void> localVarReturnType = new ParameterizedTypeReference<Void>() {};
        closeGuiRequestCreation(projectId).bodyToMono(localVarReturnType).block();
    }

    /**
     * Close GUI instance of given project-space if available.
     * Close GUI instance of given project-space if available.
     * <p><b>200</b> - OK
     * @param projectId if project-space the GUI instance is connected to.
     * @throws WebClientResponseException if an error occurs while attempting to invoke the API
     */
    public ResponseEntity<Void> closeGuiWithHttpInfo(String projectId) throws WebClientResponseException {
        ParameterizedTypeReference<Void> localVarReturnType = new ParameterizedTypeReference<Void>() {};
        return closeGuiRequestCreation(projectId).toEntity(localVarReturnType).block();
    }

    /**
     * Close GUI instance of given project-space if available.
     * Close GUI instance of given project-space if available.
     * <p><b>200</b> - OK
     * @param projectId if project-space the GUI instance is connected to.
     * @return ResponseSpec
     * @throws WebClientResponseException if an error occurs while attempting to invoke the API
     */
    public ResponseSpec closeGuiWithResponseSpec(String projectId) throws WebClientResponseException {
        return closeGuiRequestCreation(projectId);
    }
    /**
     * Open GUI instance on specified project-space and bring the GUI window to foreground.
     * Open GUI instance on specified project-space and bring the GUI window to foreground.
     * <p><b>201</b> - Created
     * @param projectId of project-space the GUI instance will connect to.
     * @param readOnly open in read-only mode.
     * @param guiParameters The guiParameters parameter
     * @throws WebClientResponseException if an error occurs while attempting to invoke the API
     */
    private ResponseSpec openGuiRequestCreation(String projectId, Boolean readOnly, GuiParameters guiParameters) throws WebClientResponseException {
        Object postBody = guiParameters;
        // verify the required parameter 'projectId' is set
        if (projectId == null) {
            throw new WebClientResponseException("Missing the required parameter 'projectId' when calling openGui", HttpStatus.BAD_REQUEST.value(), HttpStatus.BAD_REQUEST.getReasonPhrase(), null, null, null);
        }
        // create path and map variables
        final Map<String, Object> pathParams = new HashMap<String, Object>();

        pathParams.put("projectId", projectId);

        final MultiValueMap<String, String> queryParams = new LinkedMultiValueMap<String, String>();
        final HttpHeaders headerParams = new HttpHeaders();
        final MultiValueMap<String, String> cookieParams = new LinkedMultiValueMap<String, String>();
        final MultiValueMap<String, Object> formParams = new LinkedMultiValueMap<String, Object>();

        queryParams.putAll(apiClient.parameterToMultiValueMap(null, "readOnly", readOnly));

        final String[] localVarAccepts = { };
        final List<MediaType> localVarAccept = apiClient.selectHeaderAccept(localVarAccepts);
        final String[] localVarContentTypes = { 
            "application/json"
        };
        final MediaType localVarContentType = apiClient.selectHeaderContentType(localVarContentTypes);

        String[] localVarAuthNames = new String[] {  };

        ParameterizedTypeReference<Void> localVarReturnType = new ParameterizedTypeReference<Void>() {};
        return apiClient.invokeAPI("/api/projects/{projectId}/gui", HttpMethod.POST, pathParams, queryParams, postBody, headerParams, cookieParams, formParams, localVarAccept, localVarContentType, localVarAuthNames, localVarReturnType);
    }

    /**
     * Open GUI instance on specified project-space and bring the GUI window to foreground.
     * Open GUI instance on specified project-space and bring the GUI window to foreground.
     * <p><b>201</b> - Created
     * @param projectId of project-space the GUI instance will connect to.
     * @param readOnly open in read-only mode.
     * @param guiParameters The guiParameters parameter
     * @throws WebClientResponseException if an error occurs while attempting to invoke the API
     */
    public void openGui(String projectId, Boolean readOnly, GuiParameters guiParameters) throws WebClientResponseException {
        ParameterizedTypeReference<Void> localVarReturnType = new ParameterizedTypeReference<Void>() {};
        openGuiRequestCreation(projectId, readOnly, guiParameters).bodyToMono(localVarReturnType).block();
    }

    /**
     * Open GUI instance on specified project-space and bring the GUI window to foreground.
     * Open GUI instance on specified project-space and bring the GUI window to foreground.
     * <p><b>201</b> - Created
     * @param projectId of project-space the GUI instance will connect to.
     * @param readOnly open in read-only mode.
     * @param guiParameters The guiParameters parameter
     * @throws WebClientResponseException if an error occurs while attempting to invoke the API
     */
    public ResponseEntity<Void> openGuiWithHttpInfo(String projectId, Boolean readOnly, GuiParameters guiParameters) throws WebClientResponseException {
        ParameterizedTypeReference<Void> localVarReturnType = new ParameterizedTypeReference<Void>() {};
        return openGuiRequestCreation(projectId, readOnly, guiParameters).toEntity(localVarReturnType).block();
    }

    /**
     * Open GUI instance on specified project-space and bring the GUI window to foreground.
     * Open GUI instance on specified project-space and bring the GUI window to foreground.
     * <p><b>201</b> - Created
     * @param projectId of project-space the GUI instance will connect to.
     * @param readOnly open in read-only mode.
     * @param guiParameters The guiParameters parameter
     * @return ResponseSpec
     * @throws WebClientResponseException if an error occurs while attempting to invoke the API
     */
    public ResponseSpec openGuiWithResponseSpec(String projectId, Boolean readOnly, GuiParameters guiParameters) throws WebClientResponseException {
        return openGuiRequestCreation(projectId, readOnly, guiParameters);
    }
}
