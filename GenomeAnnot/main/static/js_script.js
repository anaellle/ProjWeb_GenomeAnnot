//--------------enabled tooltips
window.addEventListener('load', function() {
    const tooltipTriggerList = document.querySelectorAll('[data-bs-toggle="tooltip"]');
    const tooltipList = [...tooltipTriggerList].map(tooltipTriggerEl => new bootstrap.Tooltip(tooltipTriggerEl));
});


// --------------- Scroll button

$(document).ready(function() {
    var sizeFooter = $('footer').height();
    var sizeButton=$('#scrollTop').height();
    $('#scrollTop').css("bottom", (sizeButton+sizeFooter+10)+"px");

});


function topFunction() {
    document.body.scrollTop = 0; // For Safari
    document.documentElement.scrollTop = 0; // For Chrome, Firefox, IE and Opera
  } 


window.onscroll = function() {scrollFunction()};
function scrollFunction() {
    // When the user scrolls down 20px from the top of the document, show the button
  if (document.body.scrollTop > 50 || document.documentElement.scrollTop > 50) {
    $("#scrollTop").show();
  } else {
    $("#scrollTop").hide();
  }
}








