package de.unijena.bioinf.ms.rest.client.canopus;

import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;
import de.unijena.bioinf.fingerid.predictor_types.PredictorType;
import de.unijena.bioinf.ms.rest.client.AbstractClient;
import de.unijena.bioinf.ms.rest.model.JobUpdate;
import de.unijena.bioinf.ms.rest.model.canopus.CanopusData;
import de.unijena.bioinf.ms.rest.model.canopus.CanopusJobInput;
import de.unijena.bioinf.ms.rest.model.canopus.CanopusJobOutput;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.client.methods.HttpPost;
import org.apache.http.entity.ContentType;
import org.apache.http.entity.StringEntity;
import org.apache.http.impl.client.CloseableHttpClient;

import java.io.IOException;
import java.net.URI;
import java.nio.charset.StandardCharsets;

public class CanopusClient extends AbstractClient {

    public CanopusClient(URI serverUrl) {
        super(serverUrl);
    }

    public CanopusData getCanopusData(PredictorType predictorType, CloseableHttpClient client) throws IOException {
        return execute(client,
                () -> new HttpGet(buildVersionSpecificWebapiURI("/canopus/data.csv")
                        .setParameter("predictor", predictorType.toBitsAsString())
                        .build()),
                CanopusData::read
        );
    }

    public JobUpdate<CanopusJobOutput> postJobs(final CanopusJobInput input, CloseableHttpClient client) throws IOException {
        return executeFromJson(client,
                () -> {
                    final HttpPost post = new HttpPost(buildVersionSpecificWebapiURI("/canopus/" + CID + "/jobs").build());
                    String v = new ObjectMapper().writeValueAsString(input);
                    post.setEntity(new StringEntity(v, StandardCharsets.UTF_8));
                    post.setHeader("Content-Type", ContentType.APPLICATION_JSON.getMimeType());
                    return post;
                }, new TypeReference<>() {
                }
        );
    }
}


