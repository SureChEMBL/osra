package net.sf.osra;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.StringWriter;

import org.apache.commons.io.IOUtils;
import org.junit.Test;

/**
 * Sample usage of the library.
 * 
 * @author <a href="mailto:dmitry.katsubo@gmail.com">Dmitry Katsubo</a>
 */
public class OsraLibTest {

	@Test
	public void testGetVersion() {
		assertNotNull(OsraLib.getVersion());
	}

	@Test
	public void testProcessImage() throws IOException {
		final StringWriter writer = new StringWriter();

		int result = OsraLib.processImage(IOUtils.toByteArray(new FileInputStream(new File("test", "test.png"))),
				writer, "sdf", "inchi", true, true, true);

		System.out.println("OSRA completed with result:" + result + " structure:\n" + writer.toString() + "\n");
	}

	public static void main(String[] args) throws IOException {
		new OsraLibTest().testProcessImage();
	}
}
