package net.sourceforge.osra;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.StringWriter;
import java.io.Writer;

import org.apache.commons.io.IOUtils;

/**
 * JBI bridge for OSRA library.
 * 
 * @author <a href="mailto:dkatsubo@epo.org">Dmitry Katsubo</a>
 */
public class OsraLib {

	/**
	 * Process the given image with OSRA library. For more information see the corresponding <a
	 * href="https://sourceforge.net/apps/mediawiki/osra/index.php?title=Usage">CLI options</a>.
	 * 
	 * @param imageData
	 *            the image binary data
	 * @param outputStructureWriter
	 *            the writer to output the found structures in given format
	 * @param format
	 *            one of the formats, accepted by OpenBabel ("sdf", "smi", "can").
	 * @param outputConfidence
	 *            include confidence
	 * @param outputCoordinates
	 *            include box coordinates
	 * @param outputAvgBondLength
	 *            include average bond length
	 * @return 0, if the call succeeded or negative value in case of error
	 */
	public native int processImage(byte[] imageData, Writer outputStructureWriter, String format,
				boolean outputConfidence, boolean outputCoordinates, boolean outputAvgBondLength);

	public native String getVersion();

	static {
		System.loadLibrary("osra_java");
	}

	public static void main(String[] args) throws IOException {
		System.out.println("Processing file " + args[0]);

		final StringWriter writer = new StringWriter();

		int result = new OsraLib().processImage(IOUtils.toByteArray(new FileInputStream(args[0])), writer, "sdf", true,
					true, true);

		System.out.println("OSRA completed with result:" + result + " structure:'" + writer.toString() + "'");
	}
}
