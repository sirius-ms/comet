package de.unijena.bioinf.ms.persistence.storage;

import de.unijena.bioinf.ChemistryBase.utils.SimpleSerializers;
import de.unijena.bioinf.ms.persistence.model.core.statistics.FoldChange;
import de.unijena.bioinf.ms.persistence.model.core.tags.*;
import de.unijena.bioinf.storage.db.nosql.Index;
import de.unijena.bioinf.storage.db.nosql.Metadata;
import org.jetbrains.annotations.NotNull;

import java.io.IOException;

public interface StatsAndTaggingSupport {
    static Metadata buildMetadata() throws IOException {
        return buildMetadata(Metadata.build());
    }

    static Metadata buildMetadata(@NotNull Metadata sourceMetadata) throws IOException {
        return sourceMetadata
                .addSerialization(ValueDefinition.class, new ValueDefinition.Serializer(), new ValueDefinition.Deserializer())
                .addSerialization(ValueType.class, new SimpleSerializers.EnumAsNumberSerializer<>(), new SimpleSerializers.EnumAsNumberDeserializer<>(ValueType.class))
                .addSerialization(Tag.class, new Tag.Serializer(), new Tag.Deserializer())
                .addRepository(Tag.class,
                        Index.unique("taggedObjectId", "tagName"),
                        Index.nonUnique("tagName"),
                        Index.nonUnique("taggedObjectClass","tagName"))

//                //todo TAGS: add load of indexes for each of teh fields.
//                .addRepository(Tag.class,  Stream.concat(Stream.of(Index.unique("taggedObjectId", "tagName")),
//                        Arrays.stream(ValueType.values()).filter(ValueType::hasValue)
//                                .map(ValueType::getValueFieldName)
//                                .map(field -> Index.nonUnique("tagName", field, "taggedObjectId"))
//                ).toArray(Index[]::new))

                .addRepository(TagDefinition.class, Index.unique("tagName"), Index.nonUnique("tagType"))
                .addRepository(TagGroup.class, Index.unique("groupName"), Index.nonUnique("groupType"))

                .addRepository(FoldChange.CompoundFoldChange.class, Index.nonUnique("foreignId"))
                .addRepository(FoldChange.AlignedFeaturesFoldChange.class, Index.nonUnique("foreignId"))

                ;
    }
}
