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

function fillProjFilename(sfx) {
  //if (!wgAutoFillFilename[sfx]) {
  //  return;
 // }
  if (!document.getElementById) {
    return;
  }
  var projFileTouched = document.getElementById('wpProjFileTouched'+sfx);
  var touched = projFileTouched.value;
  if ( touched === "0" ) {
    var fname;
    // have to check this two ways??
    var uploadtype = document.getElementById('wpSourceType'+sfx);
    var uploadfile = document.getElementById('wpSourceTypeFile'+sfx);
    if ( (uploadtype != null && uploadtype.value == 'file')
        || (uploadfile != null && uploadfile.checked) ) {
      fname = document.getElementById('wpUploadFile'+sfx).value;
      if (fname.charAt(0) == '/')
      { var slash = fname.lastIndexOf('/');
        fname = fname.substring(slash+1,10000);
      }
      else if (fname.charAt(0) == '\\' || fname.indexOf(':\\') > -1)
      { var backslash = fname.lastIndexOf('\\');
        fname = fname.substring(backslash+1,10000);
      }
    } else {
      fname = document.getElementById('wpUploadFileURL'+sfx).value;
      var slash = fname.lastIndexOf('/');
      if (slash != -1)
        fname = fname.substring(slash+1,10000);
    }
    // Output result
    var destFile = document.getElementById('wpProjFilename'+sfx);
    if (destFile) {
      destFile.value = fname;
      //wgUploadWarningObj.checkNow(fname) ;
    }
  }
}
function replacePrefix(oldpfx, newpfx, lng) {
  if (lng.substring(0,oldpfx.length) == oldpfx) {
    return newpfx + lng.substring(oldpfx.length,10000);
  }
  return lng;
}
function fillDestPage(sfx) {
  if (!document.getElementById) {
    return;
  }
  var destPage = document.getElementById('wpDestPage'+sfx);
  if (destPage == null) {
    return;
  }
  // if nobody's edited the page name by hand, fill in the default suggestion
  var destPageTouched = document.getElementById('wpDestPageTouched'+sfx);
  var touched = destPageTouched.value;
  if ( touched === "0" ) {
    var projectname = document.getElementById('wpProjectName').value;
    var fname = '', orig_fname = '', orig_projectname = '';
    if (document.getElementById('wpProjFilename'+sfx) != null) {
      // for uploading source/project files: construct page name from filename
      fname = document.getElementById('wpProjFilename'+sfx).value;
    }
    if (document.getElementById('wpOrigFilename'+sfx) != null) {
      orig_fname = document.getElementById('wpOrigFilename'+sfx).value;
    }
    if (document.getElementById('wpOrigProjectname') != null) {
      orig_projectname = document.getElementById('wpOrigProjectname').value;
    }
    newValue = suggestPage(fname, projectname, orig_fname, orig_projectname);
    // Output result
    destPage.value = newValue;
      //wgUploadWarningObj.checkNow(fname) ;
    //} else {
      // for uploading a stored wiki page: keep page name except for
      // updating the project name if it changes
    //  var newValue = destPage.value;
    //  newValue = replacePrefix(wgLastProjectname,projectname,newValue);
    //  newValue = replacePrefix(wgOrigProjectname,projectname,newValue);
    //  newValue = replacePrefix('Image:'+wgLastProjectname,
    //      wgImageNamespace+':'+projectname,newValue);
    //  newValue = replacePrefix('File:'+wgLastProjectname,
    //      wgImageNamespace+':'+projectname,newValue);
    //  newValue = replacePrefix('Media:'+wgLastProjectname,
    //      'Media:'+projectname,newValue);
    //  newValue = replacePrefix('Image:'+wgOrigProjectname,
    //      wgImageNamespace+':'+projectname,newValue);
    //  newValue = replacePrefix('File:'+wgOrigProjectname,
    //      wgImageNamespace+':'+projectname,newValue);
    //  newValue = replacePrefix('Media:'+wgOrigProjectname,
    //      'Media:'+projectname,newValue);
      // Output result
    //  destPage.value = newValue;
      //wgUploadWarningObj.checkNow(fname) ;
    //}
  }
}

function redoDestPages() {
  if (!document.getElementById) {
    return;
  }
  var sfx = 1;
  while (document.getElementById('wpDestPage'+sfx)) {
    fillDestPage(sfx);
    ++sfx;
  }
  wgLastProjectname = document.getElementById('wpProjectName').value;
}

function doskip(sfx) {
  if (!document.getElementById) {
    return;
  }
  //document.getElementById('wpDestPage'+sfx).value = 'YO!';
  yesno = document.getElementById('wpSkip'+sfx).checked;
  document.getElementById('wpDestPage'+sfx).readOnly = yesno;
  if ((pf = document.getElementById('wpProjFilename'+sfx)) != null)
    pf.readOnly = yesno;
  if ((is = document.getElementById('wpIsSource'+sfx+'-0')) != null)
    is.readOnly = yesno;
  if ((is = document.getElementById('wpIsSource'+sfx+'-1')) != null)
    is.readOnly = yesno;
}

function uncheck(dstid) {
  if(!document.getElementById) return;
  document.getElementById(dstid).checked = false;
}

function toggleFilenameFiller(sfx) {
  if(!document.getElementById) return;
  var val = document.getElementById('wpProjFilename'+sfx).value;
  var touched = document.getElementById('wpProjFileTouched'+sfx);
  if ( touched != null )
    touched.value = (val == '' || val == ' ') ? "0" : "1";
  else
    wgAutoFillFilename[sfx] = (val == '' || val == ' ');
}

function togglePageFiller(sfx) {
  if(!document.getElementById) return;
  var val = document.getElementById('wpDestPage'+sfx).value;
  var touched = document.getElementById('wpDestPageTouched'+sfx);
  if ( touched != null )
    touched.value = (val == '' || val == ' ') ? "0" : "1";
  else
    wgAutoFillPage[sfx] = (val == '' || val == ' ');
}

