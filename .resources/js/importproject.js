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

function fillProjectname(sfx) {
  if (!wgAutoFillProjectname)
    return;
  if (!document.getElementById) {
    return;
  }
  var fname;
  // have to check this two ways??
  var uploadtype = document.getElementById('wpSourceType'+sfx);
  var uploadtypefile = document.getElementById('wpSourceTypeFile'+sfx);
  var uploadurl = document.getElementById('wpUploadFileURL'+sfx);
  var uploadfile = document.getElementById('wpUploadFile'+sfx);
  if ( (uploadtype != null && uploadtype.value == 'file')
      || (uploadtypefile != null && uploadtypefile.checked)
      || (uploadurl == null) ) {
    fname = uploadfile.value;
    if (fname.charAt(0) == '/')
    { var slash = fname.lastIndexOf('/');
      fname = fname.substring(slash+1,10000);
    }
    else if (fname.charAt(0) == '\\' || fname.indexOf(':\\') > -1)
    { var backslash = fname.lastIndexOf('\\');
      fname = fname.substring(backslash+1,10000);
    }
  } else if ( uploadurl != null ) {
    fname = uploadurl.value;
    var slash = fname.lastIndexOf('/');
    if (slash != -1)
      fname = fname.substring(slash+1,10000);
  } else {
    return;
  }
  if (wgPackageRegexp) {
    fname = fname.replace(wgPackageRegexp,'');
  }
  // Output result
  var dest = document.getElementById('wpProjectName'+sfx);
  if (dest) {
    dest.value = fname;
    //wgUploadWarningObj.checkNow(fname) ;
  }
}

function toggleProjectnameFiller(sfx) {
  if(!document.getElementById) return;
  var val = document.getElementById('wpProjectName'+sfx).value;
  wgAutoFillProjectname = (val == '' || val == ' ');
}

