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

package de.unijena.bioinf.ms.frontend.subtools.middleware;

import de.unijena.bioinf.ChemistryBase.utils.FileUtils;
import de.unijena.bioinf.ms.frontend.subtools.OutputOptions;
import de.unijena.bioinf.ms.frontend.subtools.Provide;
import de.unijena.bioinf.ms.frontend.subtools.RootOptions;
import de.unijena.bioinf.ms.frontend.subtools.StandaloneTool;
import de.unijena.bioinf.ms.frontend.workflow.Workflow;
import de.unijena.bioinf.ms.middleware.model.projects.ProjectInfo;
import de.unijena.bioinf.ms.middleware.service.gui.GuiService;
import de.unijena.bioinf.ms.middleware.service.projects.ProjectsProvider;
import de.unijena.bioinf.ms.properties.ParameterConfig;
import de.unijena.bioinf.projectspace.Instance;
import lombok.Getter;
import lombok.Setter;
import picocli.CommandLine;

import java.io.IOException;
import java.nio.file.Path;
import java.util.EnumSet;
import java.util.Optional;

@CommandLine.Command(name = "asService", aliases = {"rest", "REST"}, description = "EXPERIMENTAL/UNSTABLE: Starts SIRIUS as a background (REST) service that can be requested via a REST-API", versionProvider = Provide.Versions.class, mixinStandardHelpOptions = true)
public class MiddlewareAppOptions<I extends Instance> implements StandaloneTool<MiddlewareAppOptions<I>.Flow> {
    @Setter
    private ProjectsProvider<?> projectsProvider;
    @Setter
    private GuiService guiService;

    @CommandLine.Option(names = {"--port", "-p"}, description = "Specify the port on which the SIRIUS REST Service should run (Default: 8080).", defaultValue = "8080")
    private void setPort(int port) {
        System.getProperties().setProperty("server.port", String.valueOf(port));
    }

    @CommandLine.Option(names = {"--enable-rest-shutdown", "-s"}, description = "Allows to shut down the SIRIUS REST Service via a rest api call (/actuator/shutdown)", defaultValue = "false")
    private void setShutdown(boolean enableRestShutdown) {
        if (enableRestShutdown)
            System.getProperties().setProperty("management.endpoints.web.exposure.include", "health,shutdown");
        else
            System.getProperties().setProperty("management.endpoints.web.exposure.include", "health");

    }

    @CommandLine.Option(names = {"--gui", "-g"}, description = "Start GUI on specified project or on temporary project otherwise.")
    @Getter
    private boolean startGui;

    @Override
    public Flow makeWorkflow(RootOptions<?> rootOptions, ParameterConfig config) {
        return new Flow(rootOptions.getOutput());
    }

    @Override
    public int hashCode() {
        return super.hashCode();
    }

    public class Flow implements Workflow {
        private final OutputOptions output;

        private Flow(OutputOptions output) {
            this.output = output;
        }

        @Override
        public void run() {
            //do the project importing from the commandline
            Optional<Path> location = Optional.ofNullable(output).map(OutputOptions::getOutputProjectLocation);
            if (location.isPresent() || MiddlewareAppOptions.this.isStartGui()) {
                Path p = location.orElse(FileUtils.createTmpProjectSpaceLocation(null)); //todo should be part of project provider to be implementation independent
                try {
                    //open default project if given
                    ProjectInfo startPs = projectsProvider.createProject(
                                    p.getFileName().toString(),
                                    p.toAbsolutePath().toString(),
                                    EnumSet.noneOf(ProjectInfo.OptField.class), false);

                    if (isStartGui())
                        guiService.createGuiInstance(startPs.getProjectId());
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
        }
    }
}
