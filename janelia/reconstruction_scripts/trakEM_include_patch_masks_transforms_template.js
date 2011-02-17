importPackage( Packages.ini.trakem2.display );
importPackage( Packages.mpicbg.trakem2 );
importClass( Packages.java.util.HashMap );
importClass( Packages.java.awt.geom.AffineTransform );

var opener = new Opener();
var tiles = new HashMap();
var front = Display.getFront();
if ( front )
{
	var layerIterator = front.getLayer().getParent().getLayers().iterator();
	while ( layerIterator.hasNext() )
	{
		var tileIterator = layerIterator.next().getDisplayables( Patch ).iterator();
		while ( tileIterator.hasNext() )
		{
			var tile = tileIterator.next();
			tiles.put( tile.getFilePath(), tile );
			IJ.log( tile.getFilePath() + " : " + tile.getTitle() );
		}
	}

	var impMask;
	var tile;
	var tileCopy;
	
	// do this for untouched tiles
	//------------------------------------------------
	tile = tiles.get( "/home/saalfeld/tmp/8bit-grey.tif" );
	tile.setAffineTransform( new AffineTransform( m00, m10, m01, m11, m02, m12 ) );
	tile.updateMipmaps();
	
	//------------------------------------------------
	
	
	// do this for tiles that were split into pieces
	//------------------------------------------------
	impMask = opener.openImage( "/home/saalfeld/tmp/8bit-grey.mask0.tif" );
	tile = tiles.get( "/home/saalfeld/tmp/8bit-grey.tif" );
	tile.setAlphaMask( impMask.getProcessor() );
	tile.setAffineTransform( new AffineTransform( m00, m10, m01, m11, m02, m12 ) );
	tile.updateMipmaps();

	impMask = opener.openImage( "/home/saalfeld/tmp/8bit-grey.mask1.tif" );
	tileCopy = tile.clone();
	tile.getLayer().add( tileCopy );
	tileCopy.setAffineTransform( new AffineTransform( m00, m10, m01, m11, m02, m12 ) );
	tileCopy.setAlphaMask( impMask.getProcessor() );
	tileCopy.updateMipmaps();
	
	// ...
	
	//------------------------------------------------
}

