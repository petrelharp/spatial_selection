/* WorkingWiki extension for MediaWiki 1.13 and after
 * Copyright (C) 2010 Lee Worden <worden.lee@gmail.com>
 * http://lalashan.mcmaster.ca/theobio/projects/index.php/WorkingWiki
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 * http://www.gnu.org/copyleft/gpl.html
 */

/* Intuit what page a file should be imported into.
 * Before using this, set wgCapitalize, wgDefaults, wgImageNamespace,
 * wgImageLookup, and optionally wgOrigPageLocations, wgDestPageLocations, 
 * and wgPageRenames.
 * Any changes to this function should also be made to
 * wwfGeneratePageSuggestion() in misc.php.
 */
function suggestPage(dest_filename, dest_projectname, orig_filename, orig_projectname) {
  if (dest_filename == '')
    return '';
  // if the destination project already has a location for the file, use it
  if (typeof(wgDestPageLocations) != "undefined" 
      && wgDestPageLocations[dest_filename])
    return wgDestPageLocations[dest_filename]
        .replace(/\?P/g,dest_projectname).replace(/\?F/g,dest_filename);
  // or if the original project had a location for the file
  if (typeof(wgOrigPageLocations) != "undefined" 
      && wgOrigPageLocations[orig_filename])
  { var old_loc = wgOrigPageLocations[orig_filename]
      .replace(/\?P/g,orig_projectname).replace(/\?F/g,orig_filename);
    // see if that page is being uploaded to some other location
    if (typeof(wgPageRenames) != "undefined" && wgPageRenames[old_loc])
      return wgPageRenames[old_loc];
    // if not, translate it to the destination project and filename
    return wgOrigPageLocations[orig_filename]
      .replace(/\?P/g,dest_projectname).replace(/\?F/g,dest_filename);
  }
  // or, maybe the old project had it at a default location and that page
  // is being uploaded as something else
  if (typeof(wgPageRenames) != "undefined") {
    var leaf;
    var slash = filename.lastIndexOf('/');
    if (slash != -1)
      leaf = filename.substring(slash+1,10000);
    else
      leaf = null;
    var locations = new Array();
    for (i in wgDefaults) {
      var def_loc = wgDefaults[i].replace(/\?P/g,orig_projectname)
        .replace(/\?F/g,orig_filename);
      if (def_loc.substring(0,wgImageNamespace.length()+1) 
          == (wgImageNamespace+':'))
        def_loc = def_loc.replace(/\//g,'$');
      locations.push(def_loc);
      if (leaf)
      { var leaf_loc = wgDefaults[i].replace(/\?P/g,orig_projectname)
          .replace(/\?F/g,leaf);
        if (leaf_loc.substring(0,wgImageNamespace.length()+1) 
            == (wgImageNamespace+':'))
          leaf_loc = leaf_loc.replace(/\//g,'$');
        if (leaf_loc != def_loc)
          locations.push(leaf_loc);
      }
    }
    for (i in locations) {
      if (wgPageRenames[locations[i]])
        return wgPageRenames[locations[i]];
    }
  }
  // if none of these, construct the default suggestion.
  var dot = dest_filename.lastIndexOf('.');
  var image = false;
  if (dot != -1) {
    var ext = dest_filename.substring(dot+1,10000).toLowerCase();
    if (ext == 'png')
    { image = true; }
    else if (wgImportAsImageByDefault)
    { if (!wgTextLookup[ext])
      { image = true; }
    }
    else
    { if (wgImageLookup[ext])
      { image = true; }
    }
  }
  if (dest_projectname != '') {
    if (image)
      page = dest_projectname + '/' + dest_filename;
    else
      page = dest_projectname;
  } else {
    page = dest_filename;
  }
  if (wgCapitalize) { page = capitalize(page); }
  if (image) {
    page = wgImageNamespace + ':' + page.replace(/\//g, '$');
  }
  return page.replace(/ /g, '_');
}

function capitalize(string) {
  return string.charAt(0).toUpperCase() + string.substring(1,10000);
}
