import java.nio.file.*;
import java.io.IOException;
import java.nio.file.attribute.*;
import java.util.concurrent.TimeUnit;

def containerupdate() {

    File pathContainers = new File("${workflow.workDir}/singularity/");
    File[] directoryListing = pathContainers.listFiles();
    if (directoryListing != null) {

        try {
            for (File filename: directoryListing) {
                FileTime fileTime = Files.getLastModifiedTime(filename.toPath());
                int filestamp = fileTime.to(TimeUnit.SECONDS);
                int now = System.currentTimeSeconds();
                //Convert two weeks to seconds
                int weeks = TimeUnit.MINUTES.toSeconds(20160);
                int filecompare = now - weeks;

                if (filestamp < filecompare) {
                    filename.delete();
                } else {
                    continue;
                }
            }
        } catch (IOException e) {
            System.err.println("Cannot get the last modified time - " + e);
        }
    }

}
