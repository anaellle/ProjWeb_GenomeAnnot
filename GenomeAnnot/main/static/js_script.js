
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
// ---------------- Change filter of explore depending on the selection before search bar

function changeTypeResExplore(){
    const selectedValue = $("#res_type-select")[0].value;

    if(selectedValue=="genome"){
        $('#filter_expl_genome').show();
        $('#filter_expl_seq').hide();

        $('#searchbar').attr("placeholder", "Search Genome");

        
    }else if(selectedValue=="gene" || selectedValue=="prot" ){
        $('#filter_expl_genome').hide();
        $('#filter_expl_seq').show();
        if(selectedValue=="gene"){
            $('#searchbar').attr("placeholder", "Search Gene");
        } else if(selectedValue=="prot"){
            $('#searchbar').attr("placeholder", "Search Protein");
        }
    }

}

