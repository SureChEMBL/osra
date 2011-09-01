package net.sf.osra;

import java.io.Writer;

import org.apache.commons.logging.LogFactory;

/**
 * JNI bridge for OSRA library.
 * 
 * @author <a href="mailto:dmitry.katsubo@gmail.com">Dmitry Katsubo</a>
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
	 * @param embeddedFormat
	 *            format to be embedded into SDF ("inchi", "smi", "can").
	 * @param outputConfidence
	 *            include confidence
	 * @param outputCoordinates
	 *            include box coordinates
	 * @param outputAvgBondLength
	 *            include average bond length
	 * @return 0, if the call succeeded or negative value in case of error
	 */
	public static native int processImage(byte[] imageData, Writer outputStructureWriter, String format,
				String embeddedFormat, boolean outputConfidence, boolean outputCoordinates, boolean outputAvgBondLength);

	public static native String getVersion();

	static {
		try {
			System.loadLibrary("osra_java");
		}
		catch (UnsatisfiedLinkError e) {
			LogFactory.getLog(OsraLib.class).error(
						"Failed to load OSRA library from java.library.path ("
									+ System.getProperty("java.library.path") + "): " + e.getMessage());
			throw e;
		}
	}

	/*
	// Sample usage:
	public static void main(String[] args) throws IOException {
		final StringWriter writer = new StringWriter();

		int result = processImage(IOUtils.toByteArray(new FileInputStream(args[0])), writer, "sdf", "inchi", true,
					true, true);

		System.out.println("OSRA completed with result:" + result + " structure:\n" + writer.toString() + "\n");
	}
	*/
}
