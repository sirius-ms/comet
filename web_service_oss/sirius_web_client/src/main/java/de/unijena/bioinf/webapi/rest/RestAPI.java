

/*
 *
 *  This file is part of the SIRIUS library for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2013-2020 Kai Dührkop, Markus Fleischauer, Marcus Ludwig, Martin A. Hoffman and Sebastian Böcker,
 *  Chair of Bioinformatics, Friedrich-Schilller University.
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with SIRIUS. If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>
 */

package de.unijena.bioinf.webapi.rest;

import com.github.scribejava.core.model.OAuthResponseException;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.ChemistryBase.fp.CdkFingerprintVersion;
import de.unijena.bioinf.ChemistryBase.fp.MaskedFingerprintVersion;
import de.unijena.bioinf.ChemistryBase.fp.NPCFingerprintVersion;
import de.unijena.bioinf.ChemistryBase.fp.PredictionPerformance;
import de.unijena.bioinf.ChemistryBase.utils.IOFunctions;
import de.unijena.bioinf.auth.AuthService;
import de.unijena.bioinf.auth.LoginException;
import de.unijena.bioinf.canopus.CanopusResult;
import de.unijena.bioinf.chemdb.RESTDatabase;
import de.unijena.bioinf.confidence_score.svm.TrainedSVM;
import de.unijena.bioinf.fingerid.CanopusWebResultConverter;
import de.unijena.bioinf.fingerid.CovtreeWebResultConverter;
import de.unijena.bioinf.fingerid.FingerprintResult;
import de.unijena.bioinf.fingerid.FingerprintWebResultConverter;
import de.unijena.bioinf.fingerid.blast.BayesnetScoring;
import de.unijena.bioinf.fingerid.predictor_types.PredictorType;
import de.unijena.bioinf.ms.properties.PropertyManager;
import de.unijena.bioinf.rest.HttpErrorResponseException;
import de.unijena.bioinf.ms.rest.client.account.AccountClient;
import de.unijena.bioinf.ms.rest.client.canopus.CanopusClient;
import de.unijena.bioinf.ms.rest.client.chemdb.ChemDBClient;
import de.unijena.bioinf.ms.rest.client.chemdb.StructureSearchClient;
import de.unijena.bioinf.ms.rest.client.fingerid.FingerIdClient;
import de.unijena.bioinf.ms.rest.client.info.InfoClient;
import de.unijena.bioinf.ms.rest.client.jobs.JobsClient;
import de.unijena.bioinf.ms.rest.model.*;
import de.unijena.bioinf.ms.rest.model.canopus.CanopusCfData;
import de.unijena.bioinf.ms.rest.model.canopus.CanopusJobInput;
import de.unijena.bioinf.ms.rest.model.canopus.CanopusJobOutput;
import de.unijena.bioinf.ms.rest.model.canopus.CanopusNpcData;
import de.unijena.bioinf.ms.rest.model.covtree.CovtreeJobInput;
import de.unijena.bioinf.ms.rest.model.covtree.CovtreeJobOutput;
import de.unijena.bioinf.ms.rest.model.fingerid.FingerIdData;
import de.unijena.bioinf.ms.rest.model.fingerid.FingerprintJobInput;
import de.unijena.bioinf.ms.rest.model.fingerid.FingerprintJobOutput;
import de.unijena.bioinf.ms.rest.model.fingerid.TrainingData;
import de.unijena.bioinf.ms.rest.model.info.Term;
import de.unijena.bioinf.ms.rest.model.info.VersionsInfo;
import de.unijena.bioinf.ms.rest.model.license.Subscription;
import de.unijena.bioinf.ms.rest.model.license.SubscriptionConsumables;
import de.unijena.bioinf.ms.rest.model.worker.WorkerList;
import de.unijena.bioinf.ms.webapi.WebJJob;
import de.unijena.bioinf.rest.ConnectionError;
import de.unijena.bioinf.rest.ProxyManager;
import de.unijena.bioinf.storage.blob.BlobStorage;
import de.unijena.bioinf.webapi.AbstractWebAPI;
import de.unijena.bioinf.webapi.Tokens;
import org.apache.http.client.HttpClient;
import org.apache.http.client.config.RequestConfig;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.client.methods.HttpUriRequest;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.annotation.concurrent.ThreadSafe;
import java.io.IOException;
import java.net.URI;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.function.Supplier;

/**
 * Frontend WebAPI class, that represents the client to our backend rest api
 */

@ThreadSafe
public final class RestAPI extends AbstractWebAPI<RESTDatabase> {
    private static final Logger LOG = LoggerFactory.getLogger(RestAPI.class);

    private final WebJobWatcher jobWatcher = new WebJobWatcher(this);

    private final AccountClient accountClient;

    private final InfoClient serverInfoClient;
    private final JobsClient jobsClient;
    private final StructureSearchClient chemDBClient;
    private final FingerIdClient fingerprintClient;
    private final CanopusClient canopusClient;

    private static final String JOB_WATCHER_CLIENT_ID = "JOB_WATCHER";
    private static final String JOB_SUBMITTER_CLIENT_ID = "JOB_SUBMITTER";
    private static final String CONNECTION_CHECK_CLIENT_ID = "CONNECTION_CHECK";

    private Subscription activeSubscription;


    public RestAPI(@Nullable AuthService authService, @NotNull AccountClient accountClient, @NotNull InfoClient infoClient, JobsClient jobsClient, @NotNull ChemDBClient chemDBClient, @NotNull FingerIdClient fingerIdClient, @NotNull CanopusClient canopusClient) {
        super(authService);
        this.accountClient = accountClient;
        this.serverInfoClient = infoClient;
        this.jobsClient = jobsClient;
        this.chemDBClient = chemDBClient;
        this.fingerprintClient = fingerIdClient;
        this.canopusClient = canopusClient;
    }


    public RestAPI(@NotNull AuthService authService, @Nullable Subscription activeSubscription) {
        super(authService);
        IOFunctions.IOConsumer<HttpUriRequest> subsDeco = (req) -> {
            if (this.activeSubscription != null)
                req.addHeader("SUBSCRIPTION", this.activeSubscription.getSid());
        };

        this.accountClient = new AccountClient(
                URI.create(PropertyManager.getProperty("de.unijena.bioinf.sirius.web.licenseServer")),
                PropertyManager.getProperty("de.unijena.bioinf.sirius.web.licenseServer.version"),
                authService, authService, subsDeco);
        this.serverInfoClient = new InfoClient(null, authService, subsDeco);
        this.jobsClient = new JobsClient(null, authService, subsDeco);
        this.chemDBClient = new ChemDBClient(null, authService, subsDeco);
        this.fingerprintClient = new FingerIdClient(null, authService, subsDeco);
        this.canopusClient = new CanopusClient(null, authService, subsDeco);

        if (activeSubscription != null)
            changeActiveSubscription(activeSubscription);
    }


    public RestAPI(@NotNull AuthService authService, @NotNull AuthService.Token token) {
        this(authService, Tokens.getActiveSubscription(token));
    }

    public URI getSignUpURL() {
        return getAuthService().signUpURL(accountClient.getSignUpRedirectURL());
    }

    @Override
    public void changeActiveSubscription(@Nullable Subscription activeSubscription) {
        this.activeSubscription = activeSubscription;
        changeHost(this.activeSubscription != null ? () -> URI.create(this.activeSubscription.getServiceUrl()) : () -> null);
    }

    public Subscription getActiveSubscription() {
        return activeSubscription;
    }

    @Override
    public void changeHost(Supplier<URI> hostSupplier) {
        this.serverInfoClient.setServerUrl(hostSupplier);
        this.jobsClient.setServerUrl(hostSupplier);
        this.chemDBClient.setServerUrl(hostSupplier);
        this.fingerprintClient.setServerUrl(hostSupplier);
        this.canopusClient.setServerUrl(hostSupplier);
    }

    @Override
    public boolean deleteAccount() {
        return ProxyManager.doWithClient(accountClient::deleteAccount);
    }


    @Override
    public void shutdown() throws IOException {
        jobWatcher.shutdown();
        super.shutdown();
    }

    @Override
    public void acceptTermsAndRefreshToken() throws LoginException {
        if (ProxyManager.doWithClient(accountClient::acceptTerms))
            authService.refreshIfNeeded(true);
    }


    //region ServerInfo
    @NotNull
    public VersionsInfo getVersionInfo(boolean updateInfo) throws IOException {
        return ProxyManager.applyClient((c) -> serverInfoClient.getVersionInfo(c, updateInfo));
    }

    @Override
    public String getChemDbDate() {
        try {
            return ProxyManager.applyClient(chemDBClient::getChemDbDate);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public synchronized Multimap<ConnectionError.Klass, ConnectionError> checkConnection() {
        final Multimap<ConnectionError.Klass, ConnectionError> errors = Multimaps.newSetMultimap(new HashMap<>(), LinkedHashSet::new);

        try {
            ProxyManager.consumeClient(client -> {
                checkSecuredConnection(client).ifPresent(e -> errors.put(e.getErrorKlass(), e));
                //failed
                if (!errors.isEmpty()) {
                    checkLogin().ifPresent(e -> errors.put(e.getErrorKlass(), e));
                    checkUnsecuredConnection(client).ifPresent(e -> errors.put(e.getErrorKlass(), e));

                    ProxyManager.checkInternetConnection(client).ifPresent(es -> es.forEach(e -> errors.put(e.getErrorKlass(), e)));
                }
            });
        } catch (Exception e) {
            ConnectionError c = new ConnectionError(100, "Unexpected error during connection check!", ConnectionError.Klass.UNKNOWN, e);
            LOG.error(c.getSiriusMessage(), e);
            errors.put(c.getErrorKlass(), c);
        }
        return errors;
    }

    private Optional<ConnectionError> checkLogin() {
        //6,5,4
        if (authService.needsLogin())
            return Optional.of(new ConnectionError(4, "You are not logged in. Please login with a verified user account to connect to the SIRIUS web services. " +
                    "If you do not have an account you can create one for free using your institutional email address.", ConnectionError.Klass.LOGIN));
        try {
            AuthService.Token token = authService.refreshIfNeeded();

            if (!Tokens.isUserEmailVerified(token)) {
                String email = Tokens.getUserEmail(token).orElse("N/A");
                return Optional.of(new ConnectionError(51,
                        "Your accounts (primary) email address '" + email + "' has not been verified."
                                + "Please verify this email address by clicking on the verification " +
                                "link we sent to your inbox and re-login or refresh your access_token afterwards. " +
                                "Please contact support if you have not received a verification email.",
                        ConnectionError.Klass.LICENSE));
            }

            @NotNull List<Subscription> subs = Tokens.getSubscriptions(token);
            if (subs.isEmpty())
                return Optional.of(new ConnectionError(52,
                        "No Subscriptions (Licenses) found for your Account. " +
                                "Are you using your correct institutional email address? " +
                                "Subscriptions are usually bound to institutional email addresses. " +
                                "If you believe that a subscription should be available, try to re-login or refresh " +
                                "your access_token. If the problem persists contact support."
                        , ConnectionError.Klass.LICENSE));

            @Nullable Subscription sub = Tokens.getActiveSubscription(subs);
            if (sub == null)
                return Optional.of(new ConnectionError(53,
                        "Could not determine an active subscription, but there are'"
                                + subs.size() + "' subscriptions available for your account. This is likely to be a bug. " +
                                "Please contact support.", ConnectionError.Klass.LICENSE));

            java.sql.Date expDate = sub.getExpirationDate();
            if (expDate != null && expDate.getTime() < System.currentTimeMillis())
                return Optional.of(new ConnectionError(54,
                        "The active subscription expired at '"
                                + sub.getExpirationDate() + "'. Please renew your subscription or choose another non expired subscription if available." +
                                "Please contact support.", ConnectionError.Klass.LICENSE));

            @NotNull List<Term> terms = Tokens.getAcceptedTerms(token);
            String pp = sub.getPp();
            String tos = sub.getTos();
            if (tos != null && terms.stream().filter(t -> t.getLink().toString().equals(tos)).findAny().isEmpty())
                return Optional.of(new ConnectionError(61, "Terms of Service (ToS) not Accepted. Please accept the ToS of active subscription.", ConnectionError.Klass.TERMS));
            if (pp != null && terms.stream().filter(t -> t.getLink().toString().equals(pp)).findAny().isEmpty())
                return Optional.of(new ConnectionError(62, "Privacy Policy (PP) not Accepted. Please accept the PP of active subscription.", ConnectionError.Klass.TERMS));

            return Optional.empty();

        } catch (LoginException e) {
            String m = "Error when requesting login token.";
            LoggerFactory.getLogger(getClass()).error(m + ": " + e.getMessage(), e);
            return Optional.of(new ConnectionError(71, m, ConnectionError.Klass.TOKEN, e));
        }
    }

    private Optional<ConnectionError> checkUnsecuredConnection(@NotNull HttpClient client) {
        try {
            serverInfoClient.execute(client, () -> new HttpGet(serverInfoClient.getBaseURI("/actuator/health").build()));
            return Optional.empty();
        } catch (HttpErrorResponseException e) {
            String message = "Could not load version info (unsecured api endpoint). Bad Response Code.";
            LoggerFactory.getLogger(getClass()).warn(message + " Cause: " + e.getMessage());
            return Optional.of(new ConnectionError(81, message, ConnectionError.Klass.APP_SERVER, e));
        } catch (Exception e) {
            return Optional.of(new ConnectionError(82, "Unexpected error when contacting unsecured api endpoint", ConnectionError.Klass.APP_SERVER, e));
        }
    }

    private Optional<ConnectionError> checkSecuredConnection(@NotNull HttpClient client) {
        try {
            serverInfoClient.execute(client, () -> {
                HttpGet get = new HttpGet(serverInfoClient.getBaseURI("/api/check").build());
                final int timeoutInSeconds = 8000;
                get.setConfig(RequestConfig.custom().setConnectTimeout(timeoutInSeconds).setSocketTimeout(timeoutInSeconds).build());
                return get;
            });
            return Optional.empty();
        } catch (HttpErrorResponseException e) {
            String message = "Could not reach secured api endpoint. Bad Response Code.";
            LoggerFactory.getLogger(getClass()).warn(message + " Cause: " + e.getMessage());

            String[] splitMsg = e.getReasonPhrase().split(SecurityService.ERROR_CODE_SEPARATOR);
            if (splitMsg.length > 1 && splitMsg[1].equals(SecurityService.TERMS_MISSING))
                return Optional.of(new ConnectionError(63, "Server detected that that Terms and Conditions and/or Privacy policy of the used subscription has not been accepted.", ConnectionError.Klass.TERMS));

            return Optional.of(new ConnectionError(91, message, ConnectionError.Klass.APP_SERVER, e));
        } catch (OAuthResponseException e) {
            String m = "Error when contacting login Server during application server connection.";
            LoggerFactory.getLogger(getClass()).error(m + ": " + e.getMessage(), e);
            return Optional.of(new ConnectionError(72, m, ConnectionError.Klass.TOKEN, e));
        } catch (Exception e) {
            return Optional.of(new ConnectionError(92, "Unexpected error when contacting secured api endpoint", ConnectionError.Klass.APP_SERVER, e));
        }
    }

    public WorkerList getWorkerInfo() throws IOException {
        return ProxyManager.applyClient(serverInfoClient::getWorkerInfo);
    }
    //endregion

    //region Jobs
    public SubscriptionConsumables getConsumables(boolean byMonth) throws IOException {
        return getConsumables(new Date(System.currentTimeMillis()), byMonth);
    }

    public SubscriptionConsumables getConsumables(@NotNull Date monthAndYear, boolean byMonth) throws IOException {
        return ProxyManager.applyClient(client -> serverInfoClient.getConsumables(monthAndYear, byMonth, client));
    }

    public EnumMap<JobTable, List<JobUpdate<?>>> submitJobs(JobInputs submission) throws IOException {
        return ProxyManager.applyClient(client -> jobsClient.postJobs(submission, client), JOB_SUBMITTER_CLIENT_ID);
    }

    public List<JobUpdate<?>> getJobsByState(JobTable jobTable, List<JobState> statesToInclude) throws IOException {
        return getJobsByState(EnumSet.of(jobTable), statesToInclude).get(jobTable);
    }

    public EnumMap<JobTable, List<JobUpdate<?>>> getJobsByState(Collection<JobTable> jobTablesToCheck, List<JobState> statesToInclude) throws IOException {
        return ProxyManager.applyClient(client -> jobsClient.getJobsByStates(jobTablesToCheck, statesToInclude, client), JOB_WATCHER_CLIENT_ID);
    }

    public void deleteJobs(Collection<JobId> jobsToDelete, Map<JobId, Integer> countingHashes) throws IOException {
        ProxyManager.consumeClient(client -> jobsClient.deleteJobs(jobsToDelete, countingHashes, client), JOB_WATCHER_CLIENT_ID);
    }

    public void resetJobs(Collection<JobId> jobsToDelete) throws IOException {
        ProxyManager.consumeClient(client -> jobsClient.resetJobs(jobsToDelete, client), JOB_WATCHER_CLIENT_ID);
    }

    public void deleteClientAndJobs() throws IOException {
        ProxyManager.consumeClient(jobsClient::deleteAllJobs);
    }
    //endregion

    //region ChemDB
    public void consumeStructureDB(long filter, @Nullable BlobStorage cacheDir, IOFunctions.IOConsumer<RESTDatabase> doWithClient) throws IOException {
        ProxyManager.consumeClient(client -> {
            try (RESTDatabase restDB = new RESTDatabase(cacheDir, filter, chemDBClient, client)) {
                doWithClient.accept(restDB);
            }
        });
    }

    public <T> T applyStructureDB(long filter, @Nullable BlobStorage cacheDir, IOFunctions.IOFunction<RESTDatabase, T> doWithClient) throws IOException {
        return ProxyManager.applyClient(client -> {
            try (RESTDatabase restDB = new RESTDatabase(cacheDir, filter, chemDBClient, client)) {
                return doWithClient.apply(restDB);
            }
        });
    }
    //endregion

    //region Canopus
    @Override
    public WebJJob<CanopusJobInput, ?, CanopusResult, ?> submitCanopusJob(CanopusJobInput input, @Nullable Integer countingHash) throws IOException {
        final MaskedFingerprintVersion version = getClassifierMaskedFingerprintVersion(input.predictor.toCharge());
        try {
            final WebJobWatcher.SubmissionWaiterJJob<CanopusJobInput, CanopusJobOutput, CanopusResult> callback =
                    jobWatcher.submitAndWatchJob(input, JobTable.JOBS_CANOPUS, (i, id) ->
                            new RestWebJJob<>(id, input, new CanopusWebResultConverter(version,
                                    MaskedFingerprintVersion.allowAll(NPCFingerprintVersion.get()))));


            RestWebJJob<CanopusJobInput, CanopusJobOutput, CanopusResult> job = callback.awaitResult();
            job.setCountingHash(countingHash);
            return job;
        } catch (ExecutionException e) {
            throw new IOException(e);
        }
    }

    @Override
    protected CanopusCfData getCanopusCfDataUncached(@NotNull PredictorType predictorType) throws IOException {
        return ProxyManager.applyClient(client -> canopusClient.getCfData(predictorType, client));
    }

    @Override
    protected CanopusNpcData getCanopusNpcDataUncached(@NotNull PredictorType predictorType) throws IOException {
        return ProxyManager.applyClient(client -> canopusClient.getNpcData(predictorType, client));
    }
    //endregion

    //region CSI:FingerID
    public WebJJob<FingerprintJobInput, ?, FingerprintResult, ?> submitFingerprintJob(FingerprintJobInput input) throws IOException {
        final MaskedFingerprintVersion version = getCDKMaskedFingerprintVersion(input.experiment.getPrecursorIonType().getCharge());
        try {
            WebJobWatcher.SubmissionWaiterJJob<FingerprintJobInput, FingerprintJobOutput, FingerprintResult> callback =
                    jobWatcher.submitAndWatchJob(input, JobTable.JOBS_FINGERID, (in, id) ->
                            new RestWebJJob<>(id, in, new FingerprintWebResultConverter(version)));

            return callback.awaitResult();
        } catch (ExecutionException e) {
            throw new IOException(e);
        }
    }

    @Override
    protected FingerIdData getFingerIdDataUncached(@NotNull PredictorType predictorType) throws IOException {
        return ProxyManager.applyClient(client -> fingerprintClient.getFingerIdData(predictorType, client));
    }

    // use via predictor/scoring method
    public WebJJob<CovtreeJobInput, ?, BayesnetScoring, ?> submitCovtreeJob(@NotNull MolecularFormula formula, @NotNull PredictorType predictorType) throws IOException {
        final MaskedFingerprintVersion fpVersion = getFingerIdData(predictorType).getFingerprintVersion();
        final PredictionPerformance[] performances = getFingerIdData(predictorType).getPerformances();
        final CovtreeJobInput input = new CovtreeJobInput(formula.toString(), predictorType);
        try {
            WebJobWatcher.SubmissionWaiterJJob<CovtreeJobInput, CovtreeJobOutput, BayesnetScoring> callback =
                    jobWatcher.submitAndWatchJob(input, JobTable.JOBS_COVTREE, (i, id) ->
                            new RestWebJJob<>(id, i, new CovtreeWebResultConverter(fpVersion, performances)));
            return callback.awaitResult();
        } catch (ExecutionException e) {
            throw new IOException(e);
        }
    }


    /**
     * @param predictorType pos or neg
     * @param formula       Molecular formula for which the tree is requested (Default tree will be used if formula is null)
     * @return {@link BayesnetScoring} for the given {@link PredictorType} and {@link MolecularFormula}
     * @throws IOException if something went wrong with the web query
     */
    //uncached -> access via predictor
    public BayesnetScoring getBayesnetScoring(@NotNull PredictorType predictorType, @Nullable MolecularFormula formula) throws IOException {
        final MaskedFingerprintVersion fpVersion = getFingerIdData(predictorType).getFingerprintVersion();
        final PredictionPerformance[] performances = getFingerIdData(predictorType).getPerformances();
        return ProxyManager.applyClient(client -> fingerprintClient.getCovarianceScoring(predictorType, fpVersion, formula, performances, client));
    }


    //uncached -> access via predictor
    public Map<String, TrainedSVM> getTrainedConfidence(@NotNull PredictorType predictorType) throws IOException {
        return ProxyManager.applyClient(client -> fingerprintClient.getTrainedConfidence(predictorType, client));
    }

    //uncached -> access via predictor
    public TrainingData getTrainingStructures(PredictorType predictorType) throws IOException {
        return ProxyManager.applyClient(client -> fingerprintClient.getTrainingStructures(predictorType, client));
    }
    //endRegion

    //region FingerprintVersions

    /**
     * @return The Fingerprint version used by the rest Database --  not really needed but for sanity checks
     * @throws IOException if connection error happens
     */
    public CdkFingerprintVersion getCDKChemDBFingerprintVersion() throws IOException {
        return ProxyManager.applyClient(chemDBClient::getCDKFingerprintVersion);
    }
    //endregion
}
