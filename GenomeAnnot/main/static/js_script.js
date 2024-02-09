//--------------enabled tooltips
window.addEventListener('load', function() {
    const tooltipTriggerList = document.querySelectorAll('[data-bs-toggle="tooltip"]');
    const tooltipList = [...tooltipTriggerList].map(tooltipTriggerEl => new bootstrap.Tooltip(tooltipTriggerEl));
});

//----------------Hide automatic pagination

/* $(document).ready(function(){
    alert("cool");
    $("table > .pagination").hide();
}); */

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

//----------------- Create User / Login


function changeToLogin(){
    $('#formCreateUser').css("display","none");
    $('#user_left').removeClass('d-none');
    $('#user_right').addClass('d-none');
    $('#login_left').addClass('d-none');
    $('#login_right').removeClass('d-none');
    $('#titleLogin').html("Login");
}

function changeToCreateUser(){
    $('#formCreateUser').css("display","block");
    $('#user_left').addClass('d-none');
    $('#user_right').removeClass('d-none');
    $('#login_left').removeClass('d-none');
    $('#login_right').addClass('d-none');
    $('#titleLogin').html("Create User");
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

// change filter on new page after submit (when loading)
if (window.location.pathname.endsWith('/explore/')) {
    
    window.addEventListener('load', changeTypeResExplore);
}





//-------------------------- Copy Paste Sequence
/* 
function copyPaste(textElement) {
    alert("test")
    textElement.select();
    textElement.setSelectionRange(0, 99999); 
    navigator.clipboard.writeText(textElement.value);
    alert("Sequence copied !");
  } 

$('#copygene').on("click", function() {
    alert('ok');
    var seqgene=$('#seqgene');
    copyPaste(seqgene);
});
 */

