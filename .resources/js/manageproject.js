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

/**
 * JavaScript helper functions for Special:ManageProject page
 */
function fillAddPage() {
	if( !document.getElementById )
		return;
	// Output result
	var pageInput = document.getElementById('ww-add-sf-page');
	if ( pageInput )
  { var filename = document.getElementById('ww-add-sf-filename').value;
    var projectname = document.getElementById('ww-project').value;
		pageInput.value = suggestPage(filename, projectname,'','');
  }
}

function fillAddApfPage() {
	if( !document.getElementById )
		return;
	// Output result
	var pageInput = document.getElementById('ww-add-apf-page');
	if ( pageInput )
  { var filename = document.getElementById('ww-add-apf-filename').value;
    var projectname = document.getElementById('ww-project').value;
		pageInput.value = suggestPage(filename, projectname, '', '');
  }
}

function enableProjectOptionsSubmit() {
  if( !document.getElementById )
    return;
  var submit = document.getElementById('ww-project-options-submit');
  if (submit)
    submit.disabled=false;
}

function enableReadOnlyUpdate(n) {
  if( !document.getElementById )
    return;
  var submit = document.getElementById('ww-update-prereq-'+n);
  if (submit)
    submit.disabled=false;
}
