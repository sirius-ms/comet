/*
 *  This file is part of the SIRIUS Software for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2013-2020 Kai Dührkop, Markus Fleischauer, Marcus Ludwig, Martin A. Hoffman, Fleming Kretschmer, Marvin Meusel and Sebastian Böcker,
 *  Chair of Bioinformatics, Friedrich-Schiller University.
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Affero General Public License
 *  as published by the Free Software Foundation; either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License along with SIRIUS.  If not, see <https://www.gnu.org/licenses/agpl-3.0.txt>
 */

package de.unijena.bioinf.ms.gui.net;

import de.unijena.bioinf.auth.UserPortal;
import de.unijena.bioinf.ms.gui.SiriusGui;
import de.unijena.bioinf.ms.gui.actions.ActionUtils;
import de.unijena.bioinf.ms.gui.actions.SiriusActions;
import de.unijena.bioinf.ms.gui.utils.BooleanJlabel;
import de.unijena.bioinf.ms.gui.utils.TwoColumnPanel;
import de.unijena.bioinf.ms.gui.utils.loading.LoadablePanel;
import de.unijena.bioinf.ms.gui.webView.WebviewHTMLTextJPanel;
import io.sirius.ms.sdk.model.ConnectionCheck;
import io.sirius.ms.sdk.model.ConnectionError;
import io.sirius.ms.sdk.model.Subscription;
import de.unijena.bioinf.ms.properties.PropertyManager;
import org.jdesktop.swingx.JXTitledSeparator;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import javax.swing.*;
import java.awt.*;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import static de.unijena.bioinf.ms.gui.net.ConnectionChecks.toLinks;
import static io.sirius.ms.sdk.model.ConnectionError.ErrorKlassEnum;
import static io.sirius.ms.sdk.model.ConnectionError.ErrorKlassEnum.*;

public class ConnectionCheckPanel extends LoadablePanel implements PropertyChangeListener {
    public static final String READ_MORE_LICENSING = "<br> <a href=\""
            + PropertyManager.getProperty("de.unijena.bioinf.sirius.docu.url") + "/account-and-license/\">"
            + "Read more</a> about account and licensing.";

    private static String addLinkToLoginError(@Nullable String errorMessage) {
        if (errorMessage == null)
            return errorMessage;
        if (errorMessage.contains("create one"))
            return errorMessage.replace("create one", "<a href=\"" + UserPortal.signUpURL() + "\">create one</a>");
        return errorMessage;
    }


    private static final String NO_ACTIVE_SUB_MSG = "<No active subscription>";
    final BooleanJlabel internet = new BooleanJlabel();
    final BooleanJlabel authServer = new BooleanJlabel();
    final BooleanJlabel licenseServer = new BooleanJlabel();
    final BooleanJlabel fingerID = new BooleanJlabel();

    final BooleanJlabel auth = new BooleanJlabel();
    final BooleanJlabel authLicense = new BooleanJlabel();
    private final @Nullable JDialog owner;
    private final SiriusGui gui;

    JLabel authLabel = new JLabel("Authenticated ?  ");
    JLabel activeSubLabel = new JLabel("Subscription active: ?  ");
    JLabel webServiceConnectionLabel = new JLabel("Connection to SIRIUS web service ?  ");

    private JPanel resultPanel = null;
    private JXTitledSeparator subInfoSeparator = null;
    private JLabel subInfoLabel = null;

    private final TwoColumnPanel main;
    private boolean noLoginButtons;

    public ConnectionCheckPanel(@NotNull SiriusGui gui, @Nullable ConnectionCheck check) {
        this(gui, check, false);
    }
    public ConnectionCheckPanel(@NotNull SiriusGui gui, @Nullable ConnectionCheck check, boolean noLoginButtons) {
        this(null, gui, check, noLoginButtons);
    }
    public ConnectionCheckPanel(@Nullable JDialog owner, @NotNull SiriusGui gui, @Nullable ConnectionCheck check) {
        this(owner, gui, check, false);
    }
    public ConnectionCheckPanel(@Nullable JDialog owner, @NotNull SiriusGui gui, @Nullable ConnectionCheck check, boolean noLoginButtons) {
        super();
        this.owner = owner;
        this.gui = gui;
        this.noLoginButtons = noLoginButtons;
        this.main = new TwoColumnPanel(GridBagConstraints.WEST, GridBagConstraints.EAST);
        this.main.setOpaque(false);
        setContentPanel(main);
        setOpaque(false);

        main.add(new JXTitledSeparator("Connection check"), 15, false);
        main.add(new JLabel("Connection to the Internet (" + PropertyManager.getProperty("de.unijena.bioinf.sirius.web.external") + ")  "), internet, 5, false);
        main.add(new JLabel("Connection to Login Server (" + PropertyManager.getProperty("de.unijena.bioinf.sirius.security.authServer") + ")  "), authServer, 5, false);
        main.add(new JLabel("Connection to License Server (" + PropertyManager.getProperty("de.unijena.bioinf.sirius.web.licenseServer") + ")  "), licenseServer, 5, false);

        main.add(authLabel, auth, 5, false);
        main.add(activeSubLabel, authLicense, 5, false);
        main.add(webServiceConnectionLabel, fingerID, 5, false);

        main.addVerticalGlue();

        if (check == null){
            setLoading(true);
            runInBackgroundAndLoad(() -> refreshPanel(gui.getConnectionMonitor().checkConnection()));
        }else {
            refreshPanel(check);
            setLoading(false, true);
        }
    }

    public void refreshPanel(@NotNull ConnectionCheck check) {
        Optional<Subscription> sub = Optional.ofNullable(check.getLicenseInfo().getSubscription());
        String licensee = sub.map(Subscription::getSubscriberName).orElse("N/A");
        String description = sub.map(Subscription::getName).orElse(null);
        String userEmail = check.getLicenseInfo().getUserEmail();
        Set<ErrorKlassEnum> errorTypes = check.getErrors().stream().map(ConnectionError::getErrorKlass).collect(Collectors.toSet());

        internet.setState(!errorTypes.contains(INTERNET));
        authServer.setState(!errorTypes.contains(LOGIN_SERVER));
        licenseServer.setState(!errorTypes.contains(LICENSE_SERVER));
        auth.setState(userEmail != null && !errorTypes.contains(LOGIN));
        authLicense.setState(check.getLicenseInfo().getSubscription() != null && !errorTypes.contains(LICENSE));
        fingerID.setState(!errorTypes.contains(APP_SERVER));

        if (auth.isTrue()) {
            authLabel.setText(userEmail != null ? "Authenticated as '" + userEmail + "'  " : "Authenticated ?  ");
        } else {
            authLabel.setText("Not authenticated! (Or cannot verify token)  ");
        }

        activeSubLabel.setText("Subscription active: (" + sub.map(Subscription::getSid).orElse(NO_ACTIVE_SUB_MSG) + ")");
        webServiceConnectionLabel.setText("Connection to SIRIUS web service (" + sub.map(Subscription::getServiceUrl).orElse(NO_ACTIVE_SUB_MSG) + ")  ");

        //remove old states
        if (resultPanel != null)
            main.remove(resultPanel);
        if (subInfoSeparator != null)
            main.remove(subInfoSeparator);
        if (subInfoLabel != null)
            main.remove(subInfoLabel);

        //create new states
        resultPanel = createResultPanel(check);
        resultPanel.setOpaque(false);
        subInfoSeparator = new JXTitledSeparator("Subscription");
        subInfoSeparator.setOpaque(false);
        subInfoLabel = new JLabel("<html>Licensed to:  <b>" + licensee + "</b>"
                + (description != null ? " (" + description + ")" : "")
                + "</html>");
        subInfoLabel.setOpaque(false);

        //add new states
        main.add(resultPanel, 15, true);
        main.add(subInfoSeparator, 15, false);
        main.add(subInfoLabel, 5, false);

        main.revalidate();
        main.repaint();
    }


    private JPanel createResultPanel(@NotNull ConnectionCheck check) {
        TwoColumnPanel resultPanel = new TwoColumnPanel();
        resultPanel.setBorder(BorderFactory.createEmptyBorder());
        resultPanel.add(new JXTitledSeparator("Description"), 15, false);

        final List<ConnectionError> errors = check.getErrors();

        if (errors.isEmpty()  || errors.stream().filter(i -> !i.getErrorType().equals(ConnectionError.ErrorTypeEnum.WARNING)).findAny().isEmpty()) {
            //case 0 NO ERROR
            addHTMLTextPanel(resultPanel, styleGoodMessage("Connection to SIRIUS web services successfully established!"));
        } else {
            ConnectionError err = errors.getFirst();
            ErrorKlassEnum mainError = err.getErrorKlass();
            //check if internet error is not just internet check unreachable
            if (mainError == INTERNET){
                if (errors.size() > 1 &&
                        (errors.stream().map(ConnectionError::getErrorKlass).noneMatch(e -> e == LOGIN_SERVER)
                    || errors.stream().map(ConnectionError::getErrorKlass).noneMatch(e -> e == LICENSE_SERVER))
                ){
                    err = errors.stream().filter(e -> e.getErrorKlass() != INTERNET).sorted().findFirst().get();
                    mainError = err.getErrorKlass();
                }
            }
            switch (mainError) {
                case INTERNET :
                    addHTMLTextPanel(resultPanel,
                            styleErrorMessage(
                            "Could not connect to the Internet. Please check whether your computer is connected to the internet. " +
                            "All features depending on the SIRIUS web services won't work without internet connection." +
                            "If you use a proxy, please check the proxy settings." + addHtmlErrorText(err))
                    );
                    break;
                case LOGIN_SERVER:
                    addHTMLTextPanel(resultPanel,
                            styleErrorMessage(err.getSiriusMessage() + "<br>" +
                            "Could not reach " + PropertyManager.getProperty("de.unijena.bioinf.sirius.security.authServer") + ". <br>" +
                            "Either our Authentication/Login server is temporary not available " +
                            "or it cannot be reached because of your network configuration.<br>" + addHtmlErrorText(err))
                    );
                    break;
                case LICENSE_SERVER:
                    addHTMLTextPanel(resultPanel,
                            styleErrorMessage(err.getSiriusMessage() + "<br>" +
                            "Could not reach " + PropertyManager.getProperty("de.unijena.bioinf.sirius.web.licenseServer") + ". <br>" +
                            "Either our Licensing server is temporary not available, " +
                            "or it cannot be reached because of your network configuration.<br>" + addHtmlErrorText(err))
                    );
                    break;
                case TOKEN:
                    addHTMLTextPanel(resultPanel,
                            styleErrorMessage(err.getSiriusMessage() +
                            "<br>Unexpected error when refreshing/validating your access_token. <br> <b>Please try to re-login:</b>"
                            + addHtmlErrorText(err))
                    );
                    if (!noLoginButtons)
                        resultPanel.add(new JButton(ActionUtils.deriveFrom(
                            evt -> Optional.ofNullable(owner).ifPresent(JDialog::dispose),
                            SiriusActions.SIGN_OUT.getInstance(gui, true))));
                case LOGIN:
                    addHTMLTextPanel(resultPanel,
                            styleWarningMessage(addLinkToLoginError(err.getSiriusMessage()) + READ_MORE_LICENSING + addHtmlErrorText(err))
                    );
                    if (!noLoginButtons)
                        resultPanel.add(new JButton(ActionUtils.deriveFrom(
                            evt -> Optional.ofNullable(owner).ifPresent(JDialog::dispose),
                            SiriusActions.SIGN_IN.getInstance(gui, true))));
                    break;
                case LICENSE:
                    addHTMLTextPanel(resultPanel,
                            styleWarningMessage(err.getSiriusMessage() + READ_MORE_LICENSING + addHtmlErrorText(err))
                    );
                    break;
                case TERMS:
                    addHTMLTextPanel(resultPanel,
                            styleWarningMessage(err.getSiriusMessage() + "<br><b>" + toLinks(check.getLicenseInfo().getTerms()) +
                            " for the selected Webservice have not been accepted.</b><br> Click Accept to get Access:")
                    );
                    resultPanel.add(new JButton(ActionUtils.deriveFrom(
                            evt -> Optional.ofNullable(owner).ifPresent(JDialog::dispose),
                            SiriusActions.ACCEPT_TERMS.getInstance(gui, true))));
                    break;
                case APP_SERVER:
                    addHTMLTextPanel(resultPanel,
                            styleErrorMessage(err.getSiriusMessage() + "<br>" +
                            "Could not connect to the SIRIUS web service (REST API). <br><br>" + "<b>Possible reasons:</b><br>" +
                            "<ol>" +
                            "<li>The SIRIUS web service is temporary not available.</li>" +
                            "<li>The service cannot be reached because of your network configuration.</li>" +
                            "<li>Our Service is no longer available for your current SIRIUS version. <br>" +
                            "Please <a href=https://github.com/sirius-ms/sirius/releases/latest>download</a> " +
                            "the latest version of SIRIUS</li>" +
                            "</ol>" + addHtmlErrorText(err))
                    );


                    break;
                default:
                    addHTMLTextPanel(resultPanel,
                            styleErrorMessage(err.getSiriusMessage() + "<br><b>An Unexpected Network Error has occurred! " +
                            "Please submit a bug report.</b>" + addHtmlErrorText(err))
                    );
            }
        }
        return resultPanel;
    }

    private String styleErrorMessage(String errorMessage) {
        return WebviewHTMLTextJPanel.styleErrorColor(errorMessage);
    }

    private String styleWarningMessage(String errorMessage) {
        return WebviewHTMLTextJPanel.styleWarningColor(errorMessage);
    }

    private String styleGoodMessage(String errorMessage) {
        return WebviewHTMLTextJPanel.styleGoodColor(errorMessage);
    }

    private static String addHtmlErrorText(ConnectionError e){
        if (e == null || e.getServerResponseErrorCode() == null && e.getServerResponseErrorMessage() == null)
            return "";

        StringBuilder detail = new StringBuilder("<br><br> <b>Details</b><br>");

        if (e.getServerResponseErrorCode() != null)
            detail.append("Response Code: ").append(e.getServerResponseErrorCode()).append("<br>");
        if (e.getServerResponseErrorMessage() != null)
            detail.append("Response Message: ").append(e.getServerResponseErrorMessage());
        return  detail.toString();
    }
    public WebviewHTMLTextJPanel addHTMLTextPanel(@NotNull JPanel resultPanel, @NotNull String text) {
        return addHTMLTextPanel(resultPanel, text, 100);
    }

    public WebviewHTMLTextJPanel addHTMLTextPanel(@NotNull JPanel resultPanel, @NotNull String text, int height) {
        WebviewHTMLTextJPanel htmlPanel = new WebviewHTMLTextJPanel(text, getBackground());
        htmlPanel.setPreferredSize(new Dimension(main.getPreferredSize().width, height));
        resultPanel.add(htmlPanel);
        htmlPanel.load();
        return htmlPanel;
    }

    @Override
    public void propertyChange(PropertyChangeEvent evt) {
        if (evt instanceof ConnectionMonitor.ConnectionStateEvent cEvt)
            refreshPanel(cEvt.getConnectionCheck());
    }
}
